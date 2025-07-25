#' PWD-ReR: Principal Component Weighted Distance Re-Randomization
#'
#' Implements the PWD-ReR method which uses principal component analysis to weight
#' covariates by importance and applies re-randomization with distance-based balance criteria.
#'
#' @param X Covariate matrix (n x p)
#' @param pa Acceptance probability threshold for balance criterion (default: 0.2)
#' @param var_explained Proportion of variance to be explained by principal components or 'Kaiser' for Kaiser criterion (default: 0.7)
#' @param max_iterations Maximum number of re-randomization attempts (default: 1000)
#' @param seed Random seed for reproducibility (default: 2020)
#' @param balance_threshold_sims Number of simulations to estimate balance threshold (default: 500)
#'
#' @return A list containing:
#' \item{allocation}{Treatment allocation vector (0/1)}
#' \item{iterations}{Number of iterations until acceptance}
#' \item{accepted}{Whether allocation was accepted within max_iterations}
#' \item{balance_threshold}{Calculated balance threshold}
#' \item{final_balance_score}{Final weighted balance score}
#' \item{balance_history}{History of balance scores across iterations}
#'
#' @details 
#' The PWD-ReR method:
#' 1. Performs PCA to determine variable importance weights
#' 2. Creates weighted Mahalanobis-type distance metric
#' 3. Estimates balance threshold through simulation
#' 4. Re-randomizes completely until weighted balance criterion is met
#'
#' @examples
#' # Generate example data
#' X <- generate_covariates(n = 100, p = 10, rho = 0.3, seed = 123)
#' 
#' # Run PWD-ReR
#' result <- pwd_rer(X, pa = 0.1, max_iterations = 500, seed = 123)
#' 
#' # Check results
#' print(paste("Accepted:", result$accepted))
#' print(paste("Iterations:", result$iterations))
#' print(paste("Balance score:", round(result$final_balance_score, 4)))
#'
#' @export
#' @importFrom MASS mvrnorm ginv
#' @importFrom stats prcomp quantile
pwd_rer <- function(X, pa = 0.2, var_explained = 0.7, max_iterations = 1000, 
                    seed = 2020, balance_threshold_sims = 500) {
  
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  
  # Input validation
  if (pa <= 0 || pa >= 1) {
    stop("pa must be between 0 and 1")
  }
  if (max_iterations <= 0) {
    stop("max_iterations must be positive")
  }
  
  # Step 1: Calculate PCA weights for variable importance
  X_scaled <- scale(X)
  pca_result <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  
  # Determine number of components
  if (var_explained == 'Kaiser') {
    # Kaiser criterion: eigenvalues > 1
    eigenvalues <- pca_result$sdev^2
    k <- sum(eigenvalues > 1)
  } else {
    # Proportion of variance explained
    prop_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
    k <- which(prop_var >= var_explained)[1]
    if (is.na(k)) k <- p  # Use all components if threshold not reached
  }
  
  k <- max(1, min(k, p))  # Ensure k is valid
  
  # Calculate variable weights from first k principal components
  rotation_matrix <- pca_result$rotation[, 1:k, drop = FALSE]
  if (k == 1) {
    var_weights <- rotation_matrix^2
  } else {
    var_weights <- rowSums(rotation_matrix^2)
  }
  var_weights <- var_weights / sum(var_weights)
  
  # Step 2: Define weighted Mahalanobis-type distance function
  calculate_weighted_mahalanobis <- function(X_matrix, treatment_vector, weights_vector) {
    n1 <- sum(treatment_vector)
    n0 <- sum(1 - treatment_vector)
    if (n1 == 0 || n0 == 0) return(Inf)
    
    # Calculate group means
    X_treated <- X_matrix[treatment_vector == 1, , drop = FALSE]
    X_control <- X_matrix[treatment_vector == 0, , drop = FALSE]
    
    mean_treated <- colMeans(X_treated)
    mean_control <- colMeans(X_control)
    mean_diff <- mean_treated - mean_control
    
    # Create weighted covariance matrix
    X_centered <- rbind(
      sweep(X_treated, 2, mean_treated),
      sweep(X_control, 2, mean_control)
    )
    
    # Weight the covariance matrix
    cov_matrix <- (n / (n1 * n0 * (n - 1))) * t(X_matrix) %*% X_matrix
    weighted_cov <- sweep(sweep(cov_matrix, 1, sqrt(weights_vector), "*"), 2, sqrt(weights_vector), "*")
    
    # Calculate weighted Mahalanobis distance
    tryCatch({
      inv_cov <- MASS::ginv(weighted_cov)
      weighted_mean_diff <- mean_diff * sqrt(weights_vector)
      distance <- as.numeric(t(weighted_mean_diff) %*% inv_cov %*% weighted_mean_diff)
      return(distance)
    }, error = function(e) {
      # Fallback to weighted sum of squared differences
      return(sum(weights_vector * mean_diff^2))
    })
  }
  
  # Step 3: Estimate balance threshold through simulation
  threshold_simulations <- replicate(balance_threshold_sims, {
    # Generate random balanced allocation
    treatment_sim <- sample(c(rep(1, n/2), rep(0, n/2)))
    calculate_weighted_mahalanobis(X, treatment_sim, var_weights)
  })
  
  balance_threshold <- quantile(threshold_simulations, pa, na.rm = TRUE)
  
  # Step 4: Re-randomization procedure
  best_allocation <- NULL
  best_balance_score <- Inf
  balance_history <- numeric()
  accepted <- FALSE
  
  for (iteration in 1:max_iterations) {
    # Generate completely random balanced allocation
    current_allocation <- sample(c(rep(1, n/2), rep(0, n/2)))
    
    # Calculate balance score
    balance_score <- calculate_weighted_mahalanobis(X, current_allocation, var_weights)
    balance_history <- c(balance_history, balance_score)
    
    # Update best allocation
    if (balance_score < best_balance_score) {
      best_allocation <- current_allocation
      best_balance_score <- balance_score
    }
    
    # Check acceptance criterion
    if (balance_score <= balance_threshold) {
      accepted <- TRUE
      
      return(list(
        allocation = current_allocation,
        iterations = iteration,
        accepted = TRUE,
        balance_threshold = balance_threshold,
        final_balance_score = balance_score,
        balance_history = balance_history,
        k_components = k,
        var_weights = var_weights,
        method = "PWD-ReR"
      ))
    }
  }
  
  # If no allocation was accepted, return best found
  return(list(
    allocation = best_allocation,
    iterations = max_iterations,
    accepted = FALSE,
    balance_threshold = balance_threshold,
    final_balance_score = best_balance_score,
    balance_history = balance_history,
    k_components = k,
    var_weights = var_weights,
    method = "PWD-ReR"
  ))
}