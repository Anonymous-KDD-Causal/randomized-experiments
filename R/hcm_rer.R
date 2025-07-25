#' HCM-ReR: Hierarchical Clustering Matched Re-Randomization
#'
#' Implements the HCM-ReR method which extends the basic HCM method by adding
#' PCA-based variable importance weighting and re-randomization to achieve
#' excellent covariate balance.
#'
#' @param X Covariate matrix (n x p), must have even number of rows
#' @param pa Acceptance probability threshold for balance criterion (default: 0.2)
#' @param var_explained Proportion of variance to be explained by principal components or 'Kaiser' for Kaiser criterion (default: 0.7)
#' @param max_iterations Maximum number of re-randomization attempts (default: 1000)
#' @param distance_method Distance method for clustering: "euclidean", "manhattan", "maximum" (default: "euclidean")
#' @param linkage_method Linkage method for clustering: "single", "complete", "average", "ward.D2" (default: "single")
#' @param seed Random seed for reproducibility (default: 2020)
#' @param balance_threshold_sims Number of simulations to estimate balance threshold (default: 500)
#'
#' @return A list containing:
#' \item{allocation}{Treatment allocation vector (0/1)}
#' \item{iterations}{Number of iterations until acceptance}
#' \item{accepted}{Whether allocation was accepted within max_iterations}
#' \item{balance_threshold}{Calculated balance threshold}
#' \item{final_balance_score}{Final weighted balance score}
#' \item{matched_pairs}{List of matched pairs information}
#' \item{n_matched_pairs}{Number of matched pairs}
#' \item{balance_history}{History of balance scores across iterations}
#' \item{var_weights}{Variable importance weights from PCA}
#' \item{k_components}{Number of principal components used}
#'
#' @details 
#' The HCM-ReR method extends HCM by:
#' 1. Performing PCA on covariates to determine variable importance weights
#' 2. Creating weighted standardized covariates for clustering
#' 3. Using hierarchical clustering to form matched pairs
#' 4. Estimating balance threshold through simulation
#' 5. Re-randomizing treatment within pairs until weighted balance criterion is met
#'
#' @examples
#' # Generate example data
#' X <- generate_covariates(n = 100, p = 10, rho = 0.3, seed = 123)
#' 
#' # Run HCM-ReR
#' result <- hcm_rer(X, pa = 0.1, max_iterations = 500, seed = 123)
#' 
#' # Check results
#' print(paste("Accepted:", result$accepted))
#' print(paste("Iterations:", result$iterations))
#' print(paste("Balance score:", round(result$final_balance_score, 4)))
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats dist hclust prcomp quantile
hcm_rer <- function(X, pa = 0.2, var_explained = 0.7, max_iterations = 1000, 
                    distance_method = "euclidean", linkage_method = "single",
                    seed = 2020, balance_threshold_sims = 500) {
  
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  
  # Input validation
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based methods")
  }
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
  
  # Step 2: Create weighted standardized covariates for clustering
  weighted_X <- sweep(X_scaled, 2, sqrt(var_weights), "*")
  
  # Step 3: Form pairs using hierarchical clustering on weighted data
  dist_matrix <- dist(weighted_X, method = distance_method)
  hc_tree <- hclust(dist_matrix, method = linkage_method)
  
  # Extract pairs from hierarchical clustering
  pairs <- extract_pairs_from_hclust(hc_tree, n)
  
  # Step 4: Define weighted balance evaluation function
  calculate_weighted_balance <- function(X_matrix, treatment_vector, weights_vector) {
    n1 <- sum(treatment_vector)
    n0 <- sum(1 - treatment_vector)
    if (n1 == 0 || n0 == 0) return(Inf)
    
    X_std <- scale(X_matrix)
    mean_diff <- colMeans(X_std[treatment_vector == 1, , drop = FALSE]) - 
      colMeans(X_std[treatment_vector == 0, , drop = FALSE])
    sum(weights_vector * mean_diff^2)
  }
  
  # Step 5: Estimate balance threshold through simulation
  threshold_simulations <- replicate(balance_threshold_sims, {
    # Generate random balanced allocation within pairs
    treatment_sim <- rep(NA, n)
    for (pair in pairs) {
      if (length(pair) == 2) {
        treatment_sim[pair] <- sample(c(0, 1), 2)
      }
    }
    calculate_weighted_balance(X, treatment_sim, var_weights)
  })
  
  balance_threshold <- quantile(threshold_simulations, pa, na.rm = TRUE)
  
  # Step 6: Re-randomization procedure
  best_allocation <- NULL
  best_balance_score <- Inf
  balance_history <- numeric()
  accepted <- FALSE
  
  for (iteration in 1:max_iterations) {
    # Generate allocation by randomizing within pairs
    current_allocation <- rep(NA, n)
    for (pair in pairs) {
      if (length(pair) == 2) {
        current_allocation[pair] <- sample(c(0, 1), 2)
      }
    }
    
    # Skip if allocation is invalid
    if (any(is.na(current_allocation))) next
    
    # Calculate balance score
    balance_score <- calculate_weighted_balance(X, current_allocation, var_weights)
    balance_history <- c(balance_history, balance_score)
    
    # Update best allocation
    if (balance_score < best_balance_score) {
      best_allocation <- current_allocation
      best_balance_score <- balance_score
    }
    
    # Check acceptance criterion
    if (balance_score <= balance_threshold) {
      accepted <- TRUE
      
      # Create matched pairs information
      matched_pairs_info <- create_matched_pairs_info_detailed(pairs, current_allocation, X, var_weights)
      
      return(list(
        allocation = current_allocation,
        iterations = iteration,
        accepted = TRUE,
        balance_threshold = balance_threshold,
        final_balance_score = balance_score,
        matched_pairs = matched_pairs_info$matched_pairs,
        n_matched_pairs = matched_pairs_info$n_matched_pairs,
        balance_history = balance_history,
        var_weights = var_weights,
        k_components = k,
        method = "HCM-ReR"
      ))
    }
  }
  
  # If no allocation was accepted, return best found
  matched_pairs_info <- create_matched_pairs_info_detailed(pairs, best_allocation, X, var_weights)
  
  return(list(
    allocation = best_allocation,
    iterations = max_iterations,
    accepted = FALSE,
    balance_threshold = balance_threshold,
    final_balance_score = best_balance_score,
    matched_pairs = matched_pairs_info$matched_pairs,
    n_matched_pairs = matched_pairs_info$n_matched_pairs,
    balance_history = balance_history,
    var_weights = var_weights,
    k_components = k,
    method = "HCM-ReR"
  ))
}

#' Create detailed matched pairs information for HCM-ReR
#'
#' @param pairs List of pairs
#' @param allocation Treatment allocation vector
#' @param X Covariate matrix
#' @param var_weights Variable importance weights
#'
#' @return List with detailed matched pairs information
#'
#' @keywords internal
create_matched_pairs_info_detailed <- function(pairs, allocation, X, var_weights) {
  matched_pairs <- list()
  
  for (i in seq_along(pairs)) {
    pair_indices <- pairs[[i]]
    if (length(pair_indices) == 2) {
      treatments <- allocation[pair_indices]
      pair_distance <- as.numeric(dist(X[pair_indices, ]))
      
      # Calculate weighted distance
      X_weighted <- sweep(X[pair_indices, ], 2, sqrt(var_weights), "*")
      weighted_distance <- as.numeric(dist(X_weighted))
      
      matched_pairs[[i]] <- list(
        indices = pair_indices,
        treatments = treatments,
        distance = pair_distance,
        weighted_distance = weighted_distance,
        covariate_values = X[pair_indices, ],
        within_pair_diff = X[pair_indices[1], ] - X[pair_indices[2], ]
      )
    }
  }
  
  return(list(
    matched_pairs = matched_pairs,
    n_matched_pairs = length(matched_pairs)
  ))
}