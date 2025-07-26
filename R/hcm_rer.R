#' HCM_ReR: Hierarchical Clustering Matched Re-Randomization
#'
#' Implements the HCM_ReR method, which uses hierarchical clustering to form matched pairs.
#' Treatment assignments are re-randomized within these pairs to achieve excellent covariate balance 
#' based on the Mahalanobis distance. 
#'
#' @param pa Acceptance probability threshold for the Mahalanobis distance (default: 0.2)
#' @param X Covariate matrix (n x p), must have an even number of rows
#' @param n_budget Maximum number of re-randomization attempts (default: 1000)
#' @param seed Random seed for reproducibility (default: 2020)
#'
#' @return A list containing:
#' \item{allocation}{Treatment allocation vector (0/1)}
#' \item{iterations}{Number of iterations until acceptance}
#' \item{indicator}{Indicator of whether the allocation met the balance criterion (`FALSE` if accepted, `TRUE` if not accepted)}
#' \item{a}{Acceptance threshold for the Mahalanobis distance}
#' \item{matched_pairs}{List of matched pairs information (indices, distances, etc.)}
#' \item{n_matched_pairs}{Number of matched pairs created}
#' \item{n_matched_units}{Number of units successfully matched into pairs}
#' \item{balance_history}{History of Mahalanobis distances across iterations}
#'
#' @details 
#' The `HCM_ReR` method works as follows:
#' 1. **Hierarchical Clustering**: Performs hierarchical clustering on the covariate matrix (`X`) to form matched pairs.
#' 2. **Re-Randomization**: Treatment assignments are randomized within pairs.
#' 3. **Balance Evaluation**: The Mahalanobis distance is calculated for each allocation, and allocations are accepted if the distance is below a threshold determined by the chi-squared distribution.
#' 4. **Stopping Criterion**: The algorithm stops when an allocation is accepted or the maximum number of iterations (`n_budget`) is reached.
#'
#' The acceptance threshold (`a`) is derived from the chi-squared distribution with degrees of freedom equal to the number of covariates (`p`).
#'
#' @examples
#' # Generate example covariate data
#' X <- generate_covariates(n = 100, p = 5, rho = 0.3, seed = 123)
#' 
#' # Run HCM_ReR
#' result <- HCM_ReR(pa = 0.1, X = X, n_budget = 500, seed = 123)
#' 
#' # Check results
#' print(paste("Accepted:", !result$indicator))
#' print(paste("Iterations:", result$iterations))
#' print(paste("Final Mahalanobis Distance:", round(result$balance_history[result$iterations], 4)))
#'
#' @export
#' @importFrom stats dist hclust
#' 
hcm_rer <- function(pa, X, n_budget = 1000, seed = 2020) {
  set.seed(seed)
  n <- nrow(X)
  d <- ncol(X)
  a <- qchisq(pa, d)
  
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based methods")
  }
  
  # Hierarchical clustering pairing (learning from hcm_randomized)
  dist_mat <- dist(X, method = "euclidean")
  hc_tree <- hclust(dist_mat, method = "single")
  
  # Extract pairs from hierarchical clustering
  pairs <- list()
  used_indices <- logical(n)
  pair_count <- 0
  
  # Process merges to create pairs
  for (i in 1:nrow(hc_tree$merge)) {
    merge_row <- hc_tree$merge[i, ]
    
    # Get actual indices
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
  
  best_w <- NULL
  best_mdist <- Inf
  balance_history <- numeric()
  
  for (iter in 1:n_budget) {
    w <- generate_allocation_from_pairs(pairs, n)
    mdist <- maha_dist(X, w)
    
    balance_history <- c(balance_history, mdist)
    
    if (mdist < best_mdist) {
      best_w <- w
      best_mdist <- mdist
    }
    
    if (mdist <= a) {
      matched_pairs_info <- convert_pairs_to_matched_format(pairs, w, X)
      return(list(
        allocation = w, w = w, ii = iter, indicator = FALSE, a = a,
        matched_pairs = matched_pairs_info$matched_pairs,
        n_matched_pairs = matched_pairs_info$n_matched_pairs,
        n_matched_units = matched_pairs_info$n_matched_units,
        iterations = iter, balance_history = balance_history
      ))
    }
  }
  
  matched_pairs_info <- convert_pairs_to_matched_format(pairs, best_w, X)
  return(list(
    allocation = best_w, w = best_w, ii = n_budget, indicator = TRUE, a = a,
    matched_pairs = matched_pairs_info$matched_pairs,
    n_matched_pairs = matched_pairs_info$n_matched_pairs,
    n_matched_units = matched_pairs_info$n_matched_units,
    iterations = n_budget, balance_history = balance_history
  ))
}
