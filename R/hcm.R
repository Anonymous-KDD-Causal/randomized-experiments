#' Hierarchical Clustering Matching (HCM)
#' 
#' Implementation of Hierarchical Clustering Matching for experimental design.
#' This method creates matched pairs using hierarchical clustering and then
#' randomly assigns treatments within each pair.
#' 
#' @version 1.0.0

source("R/utils.R")

#' Hierarchical Clustering Matching
#'
#' Creates matched pairs using hierarchical clustering and randomly assigns
#' treatments within pairs. This method ensures perfect balance and can
#' improve precision of treatment effect estimates.
#'
#' @param X Covariate matrix (n x p). Must have even number of rows.
#' @param method Clustering linkage method (default: "single")
#' @param seed Random seed for reproducibility (default: 2020)
#' @return List containing:
#'   \item{allocation}{Treatment allocation vector (0/1)}
#'   \item{matched_pairs}{List of matched pair information}
#'   \item{n_matched_pairs}{Number of matched pairs created}
#'   \item{n_matched_units}{Number of units that were matched}
#'   \item{pair_distances}{Vector of distances for each pair}
#' @export
#' @examples
#' # Generate data
#' X <- generate_covariates(d = 5, rho = 0.3, n = 100, seed = 123)
#' 
#' # Apply HCM
#' result <- hcm(X, seed = 123)
#' 
#' # Check balance
#' maha_dist(X, result$allocation)
#' 
#' # Print matching information
#' cat("Number of pairs:", result$n_matched_pairs, "\n")
#' cat("Mean pair distance:", mean(result$pair_distances), "\n")
hcm <- function(X, linkage = "single", seed = 2020) {
  set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Check for even sample size
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based methods. Current n = ", n)
  }
  
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
  
  # Step 2: Generate single allocation from pairs
  allocation <- generate_allocation_from_pairs(pairs, n)
  
  # Calculate balance
  final_balance <- maha_dist(X, allocation)
  
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
    final_balance = final_balance,
    method = "HCM"
  ))
}
#' Evaluate HCM Results
#'
#' @param X Covariate matrix
#' @param hcm_result Result from hcm() function
#' @return List of evaluation metrics
#' @export
evaluate_hcm <- function(X, hcm_result) {
  allocation <- hcm_result$allocation
  
  # Calculate balance metrics
  maha_distance <- maha_dist(X, allocation)
  smd_values <- calc_smd(X, allocation)
  
  # Calculate pair distances
  pair_distances <- sapply(hcm_result$matched_pairs, function(p) p$distance)
  
  return(list(
    mahalanobis_distance = maha_distance,
    max_smd = max(abs(smd_values), na.rm = TRUE),
    mean_smd = mean(abs(smd_values), na.rm = TRUE),
    balance_score = mean(abs(smd_values) < 0.1, na.rm = TRUE),
    mean_pair_distance = mean(pair_distances),
    median_pair_distance = median(pair_distances)
  ))
}
