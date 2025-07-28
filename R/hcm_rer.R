#' Hierarchical Clustering Rerandomization (HCM-ReR)
#' 
#' Implementation of Hierarchical Clustering with Rerandomization.
#' This method first creates matched pairs using hierarchical clustering,
#' then applies rerandomization within the pair structure to achieve
#' better covariate balance.
#' 
#' @version 1.0.0

source("R/utils.R")

#' Hierarchical Clustering Rerandomization
#'
#' Combines hierarchical clustering matching with rerandomization.
#' First creates pairs using hierarchical clustering, then rerandomizes
#' treatment assignments within pairs until achieving acceptable balance.
#'
#' @param X Covariate matrix (n x p). Must have even number of rows.
#' @param pa Acceptance probability for rerandomization (default: 0.2)
#' @param n_budget Maximum number of rerandomization attempts (default: 1000)
#' @param linkage Clustering linkage method (default: "single")
#' @param seed Random seed for reproducibility (default: 2020)
#' @return List containing:
#'   \item{allocation}{Treatment allocation vector (0/1)}
#'   \item{matched_pairs}{List of matched pair information}
#'   \item{n_matched_pairs}{Number of matched pairs created}
#'   \item{iterations}{Number of rerandomization iterations used}
#'   \item{threshold}{Mahalanobis distance threshold used}
#'   \item{final_balance}{Final achieved balance}
#'   \item{converged}{Whether the algorithm converged within budget}
#' @export
#' @examples
#' # Generate data
#' X <- generate_covariates(d = 5, rho = 0.3, n = 100, seed = 123)
#' 
#' # Apply HCM-ReR
#' result <- hcm_rer(X, pa = 0.1, n_budget = 500, seed = 123)
#' 
#' # Check results
#' cat("Converged:", result$converged, "\n")
#' cat("Iterations:", result$iterations, "\n")
#' cat("Final balance:", result$final_balance, "\n")
hcm_rer <- function(X, pa = 0.2, n_budget = 1000, linkage = "single", seed = 2020) {
  set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Check for even sample size
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based methods. Current n = ", n)
  }
  
  # Set threshold using chi-squared distribution
  threshold <- qchisq(pa, d)
  
  # Step 1: Create pairs using hierarchical clustering
  dist_mat <- dist(X, method = "euclidean")
  hc_tree <- hclust(dist_mat, method = linkage)
  
  # Extract pairs from hierarchical clustering
  pairs <- list()
  used_indices <- logical(n)
  pair_count <- 0
  
  # Process merges to create pairs
  for (i in 1:nrow(hc_tree$merge)) {
    merge_row <- hc_tree$merge[i, ]
    
    # Get actual indices (negative values indicate original points)
    idx1 <- if (merge_row[1] < 0) abs(merge_row[1]) else NULL
    idx2 <- if (merge_row[2] < 0) abs(merge_row[2]) else NULL
    
    # If both are original points and neither is used
    if (!is.null(idx1) && !is.null(idx2) && 
        !used_indices[idx1] && !used_indices[idx2]) {
      pair_count <- pair_count + 1
      pairs[[pair_count]] <- c(idx1, idx2)
      used_indices[c(idx1, idx2)] <- TRUE
      
      if (pair_count >= n/2) break
    }
  }
  
  # Handle remaining unpaired units with greedy pairing
  remaining <- which(!used_indices)
  if (length(remaining) >= 2) {
    remaining_pairs <- create_greedy_pairs(X[remaining, , drop = FALSE])
    for (i in seq_along(remaining_pairs)) {
      pair_count <- pair_count + 1
      pairs[[pair_count]] <- remaining[remaining_pairs[[i]]]
    }
  }
  
  # Step 2: Rerandomization within pairs
  best_allocation <- NULL
  best_balance <- Inf
  balance_history <- numeric()
  
  for (iter in 1:n_budget) {
    # Generate allocation from pairs
    allocation <- generate_allocation_from_pairs(pairs, n)
    
    # Calculate balance
    current_balance <- maha_dist(X, allocation)
    balance_history <- c(balance_history, current_balance)
    
    if (current_balance < best_balance) {
      best_allocation <- allocation
      best_balance <- current_balance
    }
    
    # Check if threshold is met
    if (current_balance <= threshold) {
      # Convert pairs to matched pair format
      matched_pairs <- list()
      for (i in 1:length(pairs)) {
        pair_indices <- pairs[[i]]
        if (length(pair_indices) == 2) {
          treatments <- allocation[pair_indices]
          pair_distance <- as.numeric(dist(X[pair_indices, ]))
          
          matched_pairs[[i]] <- list(
            indices = pair_indices,
            treatments = treatments,
            distance = pair_distance
          )
        }
      }
      
      return(list(
        allocation = allocation,
        matched_pairs = matched_pairs,
        n_matched_pairs = length(matched_pairs),
        n_matched_units = length(matched_pairs) * 2,
        iterations = iter,
        threshold = threshold,
        final_balance = current_balance,
        converged = TRUE,
        balance_history = balance_history,
        method = "HCM-ReR"
      ))
    }
  }
  
  # If no allocation met threshold, return best found
  matched_pairs <- list()
  for (i in 1:length(pairs)) {
    pair_indices <- pairs[[i]]
    if (length(pair_indices) == 2) {
      treatments <- best_allocation[pair_indices]
      pair_distance <- as.numeric(dist(X[pair_indices, ]))
      
      matched_pairs[[i]] <- list(
        indices = pair_indices,
        treatments = treatments,
        distance = pair_distance
      )
    }
  }
  
  return(list(
    allocation = best_allocation,
    matched_pairs = matched_pairs,
    n_matched_pairs = length(matched_pairs),
    n_matched_units = length(matched_pairs) * 2,
    iterations = n_budget,
    threshold = threshold,
    final_balance = best_balance,
    converged = FALSE,
    balance_history = balance_history,
    method = "HCM-ReR"
  ))
}

#' Evaluate HCM-ReR Results
#'
#' @param X Covariate matrix
#' @param hcm_rer_result Result from hcm_rer() function
#' @return List of evaluation metrics
#' @export
evaluate_hcm_rer <- function(X, hcm_rer_result) {
  allocation <- hcm_rer_result$allocation
  
  # Calculate balance metrics
  maha_distance <- maha_dist(X, allocation)
  smd_values <- calc_smd(X, allocation)
  
  # Calculate pair distances
  pair_distances <- sapply(hcm_rer_result$matched_pairs, function(p) p$distance)
  
  return(list(
    mahalanobis_distance = maha_distance,
    max_smd = max(abs(smd_values), na.rm = TRUE),
    mean_smd = mean(abs(smd_values), na.rm = TRUE),
    balance_score = mean(abs(smd_values) < 0.1, na.rm = TRUE),
    mean_pair_distance = mean(pair_distances),
    median_pair_distance = median(pair_distances),
    convergence_rate = hcm_rer_result$converged,
    efficiency = hcm_rer_result$iterations / max(hcm_rer_result$balance_history, na.rm = TRUE)
  ))
}
