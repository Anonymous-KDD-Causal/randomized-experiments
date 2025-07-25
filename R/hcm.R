#' HCM: Hierarchical Clustering Matched
#'
#' Implements the basic HCM method which uses hierarchical clustering to form
#' matched pairs and randomly assigns treatments within pairs without re-randomization.
#'
#' @param X Covariate matrix (n x p), must have even number of rows
#' @param distance_method Distance method for clustering: "euclidean", "manhattan", "maximum" (default: "euclidean")
#' @param linkage_method Linkage method for clustering: "single", "complete", "average", "ward.D2" (default: "single")
#' @param seed Random seed for reproducibility (default: 2020)
#'
#' @return A list containing:
#' \item{allocation}{Treatment allocation vector (0/1)}
#' \item{matched_pairs}{List of matched pairs information}
#' \item{n_matched_pairs}{Number of matched pairs}
#' \item{method}{Method name}
#' \item{clustering_info}{Information about the clustering process}
#'
#' @details 
#' The HCM method works in three steps:
#' 1. Standardizes the covariate matrix
#' 2. Uses hierarchical clustering to form matched pairs based on covariate similarity
#' 3. Randomly assigns treatments (0/1) within each pair
#'
#' This is the base method that forms the foundation for HCM-ReR.
#'
#' @examples
#' # Generate example data
#' X <- generate_covariates(n = 100, p = 10, rho = 0.3, seed = 123)
#' 
#' # Run HCM
#' result <- hcm(X, seed = 123)
#' 
#' # Check results
#' print(paste("Number of pairs:", result$n_matched_pairs))
#' print(paste("Treatment allocation:", paste(head(result$allocation), collapse = ", ")))
#'
#' @export
#' @importFrom stats dist hclust
hcm <- function(X, distance_method = "euclidean", linkage_method = "single", seed = 2020) {
  
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  
  # Input validation
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based methods")
  }
  if (n < 2 || p < 1) {
    stop("Invalid dimensions: need at least 2 observations and 1 covariate")
  }
  
  # Step 1: Standardize covariates
  X_scaled <- scale(X)
  
  # Step 2: Perform hierarchical clustering
  dist_matrix <- dist(X_scaled, method = distance_method)
  hc_tree <- hclust(dist_matrix, method = linkage_method)
  
  # Step 3: Extract pairs from hierarchical clustering
  pairs <- extract_pairs_from_hclust(hc_tree, n)
  
  # Step 4: Randomly assign treatments within pairs
  allocation <- rep(NA, n)
  
  for (pair in pairs) {
    if (length(pair) == 2) {
      # Randomly assign 0 and 1 to the pair
      treatments <- sample(c(0, 1), 2)
      allocation[pair] <- treatments
    }
  }
  
  # Create matched pairs information
  matched_pairs_info <- create_matched_pairs_info_basic(pairs, allocation, X)
  
  return(list(
    allocation = allocation,
    matched_pairs = matched_pairs_info$matched_pairs,
    n_matched_pairs = matched_pairs_info$n_matched_pairs,
    method = "HCM",
    clustering_info = list(
      distance_method = distance_method,
      linkage_method = linkage_method,
      tree_height = max(hc_tree$height)
    )
  ))
}

#' Extract pairs from hierarchical clustering tree
#'
#' @param hc_tree Hierarchical clustering object
#' @param n Total number of observations
#'
#' @return List of pairs
#'
#' @keywords internal
extract_pairs_from_hclust <- function(hc_tree, n) {
  pairs <- list()
  used_indices <- logical(n)
  pair_count <- 0
  
  # Process merge matrix to find pairs of individual observations
  for (i in 1:nrow(hc_tree$merge)) {
    merge_row <- hc_tree$merge[i, ]
    
    # Convert negative indices to positive (leaf nodes)
    idx1 <- if (merge_row[1] < 0) abs(merge_row[1]) else NULL
    idx2 <- if (merge_row[2] < 0) abs(merge_row[2]) else NULL
    
    # If both are leaf nodes and not used, create pair
    if (!is.null(idx1) && !is.null(idx2) && 
        !used_indices[idx1] && !used_indices[idx2]) {
      pair_count <- pair_count + 1
      pairs[[pair_count]] <- c(idx1, idx2)
      used_indices[c(idx1, idx2)] <- TRUE
      
      if (pair_count >= n/2) break
    }
  }
  
  # Handle any remaining unpaired units with greedy matching
  remaining <- which(!used_indices)
  if (length(remaining) >= 2) {
    remaining_pairs <- create_greedy_pairs_simple(remaining)
    for (i in seq_along(remaining_pairs)) {
      if (length(remaining_pairs[[i]]) == 2) {
        pair_count <- pair_count + 1
        pairs[[pair_count]] <- remaining_pairs[[i]]
      }
    }
  }
  
  return(pairs)
}

#' Simple greedy pairing for remaining units
#'
#' @param indices Vector of indices to be paired
#'
#' @return List of pairs
#'
#' @keywords internal
create_greedy_pairs_simple <- function(indices) {
  pairs <- list()
  remaining <- indices
  pair_count <- 0
  
  while (length(remaining) >= 2) {
    pair_count <- pair_count + 1
    pairs[[pair_count]] <- remaining[1:2]
    remaining <- remaining[-(1:2)]
  }
  
  return(pairs)
}

#' Create matched pairs information for basic HCM
#'
#' @param pairs List of pairs
#' @param allocation Treatment allocation vector
#' @param X Covariate matrix
#'
#' @return List with matched pairs information
#'
#' @keywords internal
create_matched_pairs_info_basic <- function(pairs, allocation, X) {
  matched_pairs <- list()
  
  for (i in seq_along(pairs)) {
    pair_indices <- pairs[[i]]
    if (length(pair_indices) == 2) {
      treatments <- allocation[pair_indices]
      pair_distance <- as.numeric(dist(X[pair_indices, ]))
      
      matched_pairs[[i]] <- list(
        indices = pair_indices,
        treatments = treatments,
        distance = pair_distance,
        covariate_values = X[pair_indices, ]
      )
    }
  }
  
  return(list(
    matched_pairs = matched_pairs,
    n_matched_pairs = length(matched_pairs)
  ))
}