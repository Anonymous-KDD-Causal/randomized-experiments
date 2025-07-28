#' PCA-Weighted Distance Rerandomization (PWD-ReR)
#' 
#' Implementation of PCA-Weighted Distance Rerandomization for experimental design.
#' This method uses PCA to weight variables by their importance and applies
#' rerandomization using a weighted balance criterion.
#' 
#' @version 1.0.0

source("R/utils.R")

#' PCA-Weighted Distance Rerandomization
#'
#' Uses Principal Component Analysis to weight variables and applies
#' rerandomization with a weighted balance criterion. This method is
#' particularly effective when covariates have different importance levels.
#'
#' @param X Covariate matrix (n x p)
#' @param pa Acceptance probability for rerandomization (default: 0.2)
#' @param var_explained Proportion of variance to explain with PCA (default: 0.7)
#' @param n_budget Maximum number of rerandomization attempts (default: 1000)
#' @param n_sim_threshold Number of simulations for threshold estimation (default: 500)
#' @param seed Random seed for reproducibility (default: 2020)
#' @return List containing:
#'   \item{allocation}{Treatment allocation vector (0/1)}
#'   \item{var_weights}{PCA-derived variable weights}
#'   \item{k}{Number of principal components used}
#'   \item{iterations}{Number of rerandomization iterations used}
#'   \item{threshold}{Weighted balance threshold used}
#'   \item{final_balance}{Final achieved weighted balance}
#'   \item{converged}{Whether the algorithm converged within budget}
#' @export
#' @examples
#' # Generate data with different variable importance
#' X <- generate_covariates(d = 10, rho = 0.5, n = 100, seed = 123)
#' 
#' # Apply PWD-ReR
#' result <- pwd_rer(X, pa = 0.1, var_explained = 0.8, seed = 123)
#' 
#' # Check results
#' cat("Converged:", result$converged, "\n")
#' cat("PCs used:", result$k, "\n")
#' cat("Variable weights:", round(result$var_weights, 3), "\n")
pwd_rer <- function(X, pa = 0.2, var_explained = 0.7, n_budget = 1000, 
                    n_sim_threshold = 500, seed = 2020) {
  set.seed(seed)
  
  n <- nrow(X)
  p <- ncol(X)
  n0 <- ceiling(n/2)
  n1 <- n - n0
  
  # Step 1: Calculate PCA weights
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_svd <- svd(X_centered)
  
  # Determine number of components to use
  if (var_explained == 'Kaiser') {
    # Kaiser criterion: eigenvalues > average eigenvalue
    k <- sum(X_svd$d^2 > mean(X_svd$d^2))
  } else {
    # Proportion of variance explained
    cumsum_vars <- cumsum(X_svd$d^2)
    cumsum_vars <- cumsum_vars / cumsum_vars[length(cumsum_vars)]
    k <- sum(cumsum_vars < var_explained) + 1
  }
  
  # Calculate variable weights from PCA loadings
  rotation <- X_svd$v[, 1:k, drop = FALSE]
  if (k == 1) {
    var_weights <- rotation^2
  } else {
    var_weights <- rowSums(rotation^2)
  }
  var_weights <- var_weights / sum(var_weights)  # Normalize
  
  # Step 2: Define weighted balance criterion
  calc_weighted_balance <- function(X, W, weights) {
    n1 <- sum(W)
    n0 <- sum(1 - W)
    if (n1 == 0 || n0 == 0) return(Inf)
    
    # Standardize X for fair comparison across variables
    X_std <- scale(X)
    mean_diff <- colMeans(X_std[W == 1, , drop = FALSE]) - 
      colMeans(X_std[W == 0, , drop = FALSE])
    
    # Weighted sum of squared standardized mean differences
    weighted_balance <- sum(weights * mean_diff^2)
    return(weighted_balance)
  }
  
  # Step 3: Simulate threshold for weighted balance criterion
  sim_distances <- replicate(n_sim_threshold, {
    W_sim <- sample(c(rep(1, n1), rep(0, n0)))
    calc_weighted_balance(X, W_sim, var_weights)
  })
  threshold <- quantile(sim_distances, pa, na.rm = TRUE)
  
  # Step 4: Rerandomization
  best_balance <- Inf
  best_allocation <- NULL
  balance_history <- numeric()
  
  for (iter in 1:n_budget) {
    # Generate balanced allocation
    allocation <- sample(c(rep(1, n1), rep(0, n0)))
    
    # Calculate weighted balance
    balance_score <- calc_weighted_balance(X, allocation, var_weights)
    balance_history <- c(balance_history, balance_score)
    
    # Update best allocation
    if (balance_score < best_balance) {
      best_balance <- balance_score
      best_allocation <- allocation
    }
    
    # Check convergence
    if (balance_score <= threshold) {
      return(list(
        allocation = allocation,
        var_weights = var_weights,
        k = k,
        iterations = iter,
        threshold = threshold,
        final_balance = balance_score,
        converged = TRUE,
        balance_history = balance_history[1:iter],
        method = "PWD-ReR",
        pca_rotation = rotation
      ))
    }
  }
  
  # Return best allocation if no convergence
  return(list(
    allocation = best_allocation,
    var_weights = var_weights,
    k = k,
    iterations = n_budget,
    threshold = threshold,
    final_balance = best_balance,
    converged = FALSE,
    balance_history = balance_history,
    method = "PWD-ReR",
    pca_rotation = rotation
  ))
}

#' Evaluate PWD-ReR Results
#'
#' @param X Covariate matrix
#' @param pwd_rer_result Result from pwd_rer() function
#' @return List of evaluation metrics
#' @export
evaluate_pwd_rer <- function(X, pwd_rer_result) {
  allocation <- pwd_rer_result$allocation
  weights <- pwd_rer_result$var_weights
  
  # Standard balance metrics
  maha_distance <- maha_dist(X, allocation)
  smd_values <- calc_smd(X, allocation)
  
  # Weighted balance metrics
  X_std <- scale(X)
  mean_diff <- colMeans(X_std[allocation == 1, , drop = FALSE]) - 
    colMeans(X_std[allocation == 0, , drop = FALSE])
  weighted_balance <- sum(weights * mean_diff^2)
  
  return(list(
    mahalanobis_distance = maha_distance,
    weighted_balance = weighted_balance,
    max_smd = max(abs(smd_values), na.rm = TRUE),
    mean_smd = mean(abs(smd_values), na.rm = TRUE),
    weighted_smd = sqrt(sum(weights * smd_values^2)),
    balance_score = mean(abs(smd_values) < 0.1, na.rm = TRUE),
    convergence_rate = pwd_rer_result$converged,
    pca_components = pwd_rer_result$k,
    variable_importance = pwd_rer_result$var_weights
  ))
}

#' Plot Variable Weights from PWD-ReR
#'
#' @param pwd_rer_result Result from pwd_rer() function
#' @param var_names Optional variable names
#' @return ggplot object (if ggplot2 is available)
#' @export
plot_variable_weights <- function(pwd_rer_result, var_names = NULL) {
  weights <- pwd_rer_result$var_weights
  
  if (is.null(var_names)) {
    var_names <- paste0("X", 1:length(weights))
  }
  
  # Create simple base R plot if ggplot2 not available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    barplot(weights, names.arg = var_names, 
            main = "Variable Weights from PCA",
            ylab = "Weight", xlab = "Variable",
            las = 2)
    return(invisible(NULL))
  }
  
  # Create ggplot if available
  df <- data.frame(
    Variable = factor(var_names, levels = var_names[order(weights, decreasing = TRUE)]),
    Weight = weights
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Weight)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Variable Weights from PCA",
                  y = "Weight", x = "Variable")
}
