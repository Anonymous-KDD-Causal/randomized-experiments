#' Hierarchical Clustering Matching (HCM)
#' 
#' Implementation of Hierarchical Clustering Matching for experimental design.
#' This method creates matched pairs using hierarchical clustering and then
#' randomly assigns treatments within each pair.
#' 
#' @author Your Name
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
hcm <- function(X, method = "single", seed = 2020) {
  set.seed(seed)
  
  n <- nrow(X)
  
  # Check for even sample size
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based matching. Current n = ", n)
  }
  
  original_indices <- 1:n
  matched_pairs <- list()
  pair_distances <- numeric()
  pair_count <- 0
  
  working_data <- X
  working_indices <- original_indices
  
  # Iteratively find closest pairs and match them
  while (nrow(working_data) >= 2) {
    if (nrow(working_data) == 2) {
      # Last pair
      pair_count <- pair_count + 1
      treatment_assignment <- sample(c(0, 1), 2, replace = FALSE)
      pair_distance <- as.numeric(dist(working_data))
      
      matched_pairs[[pair_count]] <- list(
        indices = working_indices,
        treatments = treatment_assignment,
        distance = pair_distance
      )
      pair_distances[pair_count] <- pair_distance
      break
    }
    
    # Build hierarchical clustering tree
    dist_mat <- dist(working_data, method = "euclidean")
    hc_tree <- hclust(dist_mat, method = method)
    
    # Find the closest pair (first merge in dendrogram)
    first_merge <- hc_tree$merge[1, ]
    
    # Extract indices of closest pair
    if (first_merge[1] < 0) {
      idx1 <- abs(first_merge[1])
    } else {
      idx1 <- 1  # Fallback
    }
    
    if (first_merge[2] < 0) {
      idx2 <- abs(first_merge[2])
    } else {
      idx2 <- 2  # Fallback
    }
    
    # Safety checks
    if (idx1 > nrow(working_data)) idx1 <- 1
    if (idx2 > nrow(working_data)) idx2 <- 2
    if (idx1 == idx2) idx2 <- min(idx1 + 1, nrow(working_data))
    
    # Create matched pair
    pair_count <- pair_count + 1
    treatment_assignment <- sample(c(0, 1), 2, replace = FALSE)
    pair_distance <- as.numeric(dist(working_data[c(idx1, idx2), ]))
    
    matched_pairs[[pair_count]] <- list(
      indices = working_indices[c(idx1, idx2)],
      treatments = treatment_assignment,
      distance = pair_distance
    )
    pair_distances[pair_count] <- pair_distance
    
    # Remove matched units from working data
    remove_indices <- c(idx1, idx2)
    working_data <- working_data[-remove_indices, , drop = FALSE]
    working_indices <- working_indices[-remove_indices]
  }
  
  # Create final allocation vector
  allocation <- rep(NA, n)
  for (i in 1:pair_count) {
    pair_info <- matched_pairs[[i]]
    allocation[pair_info$indices] <- pair_info$treatments
  }
  
  # Handle any unmatched units (should not happen with even n)
  unmatched <- which(is.na(allocation))
  if (length(unmatched) > 0) {
    warning("Found unmatched units. Assigning random treatments.")
    allocation[unmatched] <- sample(c(0, 1), length(unmatched), replace = TRUE)
  }
  
  return(list(
    allocation = allocation,
    matched_pairs = matched_pairs,
    n_matched_pairs = pair_count,
    n_matched_units = pair_count * 2,
    pair_distances = pair_distances,
    method = "HCM",
    linkage_method = method
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
  
  return(list(
    mahalanobis_distance = maha_distance,
    max_smd = max(abs(smd_values), na.rm = TRUE),
    mean_smd = mean(abs(smd_values), na.rm = TRUE),
    balance_score = mean(abs(smd_values) < 0.1, na.rm = TRUE),
    mean_pair_distance = mean(hcm_result$pair_distances),
    median_pair_distance = median(hcm_result$pair_distances)
  ))
}
