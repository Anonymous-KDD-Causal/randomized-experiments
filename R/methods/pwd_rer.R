#' PCA-Weighted Distance Rerandomization (PWD-ReR)
#' 
#' Implements rerandomization with a PCA-weighted balance criterion that gives
#' higher importance to covariates contributing more to the first few principal
#' components. This method provides targeted balance improvement on the most
#' important dimensions of the covariate space.
#' 
#' @param pa Acceptance probability for rerandomization threshold
#' @param X Matrix of covariates (n x p)
#' @param var_explained Proportion of variance explained by PCA components (default: 0.7)
#' @param n_budget Maximum number of rerandomization iterations (default: 1000)
#' @param seed Random seed for reproducibility
#' @param n_sim_threshold Number of simulations for threshold estimation (default: 500)
#' 
#' @return List containing:
#'   - allocation: Treatment allocation vector (0/1)
#'   - w: Treatment allocation vector (alias for allocation)
#'   - ii: Number of iterations used
#'   - indicator: TRUE if budget exhausted, FALSE if threshold met
#'   - a: Acceptance threshold used
#'   - k: Number of PCA components used
#'   - var_weights: PCA-derived variable importance weights
#'   - iterations: Number of iterations used (alias for ii)
#'   - balance_history: Vector of balance scores across iterations
#'   - method: Method name identifier
#' 
#' @details
#' The PWD-ReR method works by:
#' 1. Performing PCA on the covariate matrix
#' 2. Calculating variable importance weights based on loadings of first k components
#' 3. Defining a weighted balance metric that emphasizes important variables
#' 4. Estimating acceptance threshold via simulation
#' 5. Performing rerandomization until threshold is met or budget exhausted
#' 
#' The weighted balance metric is defined as:
#' weighted_balance = sum(weights * (standardized_mean_differences)^2)
#' 
#' This approach is particularly effective when covariates have different
#' importance levels or when the covariate space has a natural low-dimensional
#' structure captured by the first few principal components.
#' 
#' @examples
#' # Generate example data with factor structure
#' set.seed(2020)
#' n <- 200
#' p <- 15
#' 
#' # Create data with clear principal component structure
#' factor1 <- rnorm(n)
#' factor2 <- rnorm(n)
#' factor3 <- rnorm(n)
#' 
#' X <- cbind(
#'   factor1 + 0.3*rnorm(n), factor1 + 0.3*rnorm(n), factor1 + 0.3*rnorm(n),
#'   factor2 + 0.3*rnorm(n), factor2 + 0.3*rnorm(n), factor2 + 0.3*rnorm(n),
#'   factor3 + 0.3*rnorm(n), factor3 + 0.3*rnorm(n), factor3 + 0.3*rnorm(n),
#'   matrix(rnorm(n*6), n, 6)  # noise variables
#' )
#' X <- scale(X)
#' 
#' # Run PWD-ReR
#' result <- PWD_ReR(pa = 0.05, X = X)
#' 
#' # Check variable weights
#' print("Variable importance weights:")
#' print(round(result$var_weights, 3))
#' 
#' # Compare with standard rerandomization
#' result_standard <- rerandomization_standard(pa = 0.05, X = X)
#' cat("PWD-ReR iterations:", result$iterations, "\n")
#' cat("Standard ReR iterations:", result_standard$iterations, "\n")
#' 
#' @export
PWD_ReR <- function(pa, X, var_explained = 0.7, n_budget = 1000, 
                    seed = 2020, n_sim_threshold = 500) {
  
  # Input validation
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (pa <= 0 || pa >= 1) {
    stop("pa must be between 0 and 1")
  }
  if (n_budget <= 0) {
    stop("n_budget must be positive")
  }
  if (var_explained <= 0 || var_explained > 1) {
    stop("var_explained must be between 0 and 1")
  }
  
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  n0 <- ceiling(n/2)
  n1 <- n - n0
  
  # Step 1: Perform PCA and calculate variable importance weights
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_svd <- svd(X_centered)
  
  if (var_explained == 'Kaiser') {
    # Kaiser criterion: keep components with eigenvalues > mean eigenvalue
    k <- sum(X_svd$d^2 > mean(X_svd$d^2))
  } else {
    # Cumulative variance explained criterion
    cumsum_vars <- cumsum(X_svd$d^2)
    cumsum_vars <- cumsum_vars / cumsum_vars[length(cumsum_vars)]
    k <- sum(cumsum_vars < var_explained) + 1
  }
  
  # Ensure we have at least 1 component and at most p components
  k <- max(1, min(k, p))
  
  # Calculate variable importance weights from PCA loadings
  rotation <- X_svd$v[, 1:k, drop = FALSE]
  
  if (k == 1) {
    var_weights <- rotation^2
  } else {
    var_weights <- rowSums(rotation^2)
  }
  var_weights <- var_weights / sum(var_weights)
  
  # Step 2: Define weighted balance metric
  calc_weighted_balance <- function(X, W, weights) {
    n1 <- sum(W)
    n0 <- sum(1-W)
    if (n1 == 0 || n0 == 0) return(Inf)
    
    # Standardize covariates for fair comparison across variables
    X_std <- scale(X)
    mean_diff <- colMeans(X_std[W == 1, , drop = FALSE]) - 
      colMeans(X_std[W == 0, , drop = FALSE])
    
    # Weighted sum of squared standardized mean differences
    weighted_balance <- sum(weights * mean_diff^2)
    return(weighted_balance)
  }
  
  # Step 3: Estimate acceptance threshold via simulation
  sim_distances <- replicate(n_sim_threshold, {
    W_sim <- sample(c(rep(1, n1), rep(0, n0)))
    calc_weighted_balance(X, W_sim, var_weights)
  })
  
  threshold <- quantile(sim_distances, pa, na.rm = TRUE)
  
  # Step 4: Perform rerandomization
  best_balance <- Inf
  best_W <- NULL
  balance_history <- numeric()
  
  for(iter in 1:n_budget) {
    # Generate random allocation
    W <- sample(c(rep(1, n1), rep(0, n0)))
    
    # Calculate weighted balance score
    balance_score <- calc_weighted_balance(X, W, var_weights)
    balance_history <- c(balance_history, balance_score)
    
    # Update best allocation found so far
    if(balance_score < best_balance) {
      best_balance <- balance_score
      best_W <- W
    }
    
    # Check if threshold is met
    if(balance_score <= threshold) {
      # Accept this allocation
      return(list(
        allocation = W,
        w = W,
        indicator = FALSE,
        ii = iter,
        k = k,
        a = threshold,
        var_weights = as.vector(var_weights),
        iterations = iter,
        balance_history = balance_history,
        method = "PWD-ReR"
      ))
    }
  }
  
  # Return best allocation if no allocation was accepted
  return(list(
    allocation = best_W,
    w = best_W,
    indicator = TRUE,
    ii = n_budget,
    k = k,
    a = threshold,
    var_weights = as.vector(var_weights),
    iterations = n_budget,
    balance_history = balance_history,
    method = "PWD-ReR"
  ))
}


#' Standard Rerandomization (for comparison)
#' 
#' Implements standard rerandomization using Mahalanobis distance
#' 
#' @param pa Acceptance probability
#' @param X Covariate matrix
#' @param n_budget Maximum iterations
#' @param seed Random seed
#' 
#' @export
rerandomization_standard <- function(pa, X, n_budget = 1000, seed = 2020) {
  set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  a <- qchisq(pa, d)
  
  n0 <- ceiling(n/2)
  n1 <- n - n0
  
  # Helper function to calculate Mahalanobis distance
  maha_dist <- function(X, w, lambda = 0) {
    n1 <- sum(w)
    n0 <- sum(1-w)
    n <- n0 + n1
    xbar1 <- t(X) %*% w / n1
    xbar0 <- t(X) %*% (1-w) / n0
    delta <- xbar1 - xbar0
    cov_mat <- (n/(n1*n0*(n-1))) * t(X) %*% X
    mdist <- as.numeric(t(delta) %*% MASS::ginv(cov_mat + lambda*diag(nrow(cov_mat))) %*% delta)
    return(mdist)
  }
  
  w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
  mdist <- ifelse(n < d, n-1, maha_dist(X, w))
  
  ii <- 1
  best_w <- w
  best_mdist <- mdist
  balance_history <- c(mdist)
  
  while (mdist > a) {
    w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
    mdist <- maha_dist(X, w)
    balance_history <- c(balance_history, mdist)
    
    if(best_mdist > mdist){
      best_w <- w
      best_mdist <- mdist
    }
    ii <- ii + 1
    if(ii >= n_budget){
      break
    }
  }
  
  indicator <- ii >= n_budget
  if(indicator == TRUE){
    return_w <- best_w
  } else {
    return_w <- w
  }
  
  return(list(
    allocation = return_w,
    w = return_w, 
    ii = ii, 
    indicator = indicator, 
    a = a,
    iterations = ii,
    balance_history = balance_history,
    method = "Standard-ReR"
  ))
}


#' Visualize PWD-ReR Results
#' 
#' Creates diagnostic plots for PWD-ReR including variable weights and convergence
#' 
#' @param result PWD-ReR result object
#' @param X Original covariate matrix (optional, for variable names)
#' @param show_weights Whether to show variable importance weights plot
#' @param show_convergence Whether to show convergence plot
#' 
#' @export
visualize_PWD_ReR <- function(result, X = NULL, show_weights = TRUE, show_convergence = TRUE) {
  if (!require(ggplot2)) {
    stop("ggplot2 package required for visualization")
  }
  
  plots <- list()
  
  # Variable weights plot
  if (show_weights && !is.null(result$var_weights)) {
    p <- length(result$var_weights)
    var_names <- if(!is.null(X) && !is.null(colnames(X))) {
      colnames(X)
    } else {
      paste0("X", 1:p)
    }
    
    weights_df <- data.frame(
      variable = factor(var_names, levels = var_names),
      weight = result$var_weights,
      rank = rank(-result$var_weights)
    )
    
    p1 <- ggplot(weights_df, aes(x = variable, y = weight)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "PCA-Derived Variable Importance Weights",
           subtitle = paste("Based on first", result$k, "principal components"),
           x = "Variables", y = "Weight")
    
    plots$weights <- p1
    print(p1)
  }
  
  # Convergence plot
  if (show_convergence && !is.null(result$balance_history)) {
    balance_df <- data.frame(
      iteration = 1:length(result$balance_history),
      balance_score = result$balance_history
    )
    
    p2 <- ggplot(balance_df, aes(x = iteration, y = balance_score)) +
      geom_line(color = "blue", alpha = 0.7) +
      geom_point(size = 0.5, alpha = 0.5) +
      theme_minimal() +
      labs(title = "PWD-ReR Convergence",
           subtitle = paste("Converged in", result$iterations, "iterations"),
           x = "Iteration", y = "Weighted Balance Score")
    
    if (!is.null(result$a)) {
      p2 <- p2 + geom_hline(yintercept = result$a, color = "red", linetype = "dashed",
                            alpha = 0.7) +
        annotate("text", x = length(result$balance_history) * 0.8, y = result$a,
                 label = paste("Threshold =", round(result$a, 4)), 
                 vjust = -0.5, color = "red")
    }
    
    plots$convergence <- p2
    print(p2)
  }
  
  # Print summary
  cat("PWD-ReR Results Summary:\n")
  cat(paste("Total iterations:", result$iterations, "\n"))
  cat(paste("Budget exhausted:", result$indicator, "\n"))
  cat(paste("PCA components used:", result$k, "\n"))
  if (!is.null(result$balance_history)) {
    cat(paste("Final balance score:", round(tail(result$balance_history, 1), 6), "\n"))
    cat(paste("Acceptance threshold:", round(result$a, 6), "\n"))
    cat(paste("Improvement over iteration 1:", 
              round((result$balance_history[1] - tail(result$balance_history, 1)) / 
                      result$balance_history[1] * 100, 2), "%\n"))
  }
  
  invisible(plots)
}