#' Hierarchical Clustering Matching (HCM)
#' 
#' Implements hierarchical clustering-based matching for treatment allocation.
#' This method creates pairs of similar units using hierarchical clustering
#' and then randomly assigns treatments within each pair.
#' 
#' @param X Matrix of covariates (n x p)
#' @param seed Random seed for reproducibility
#' @param method Linkage method for hierarchical clustering (default: "single")
#' 
#' @return List containing:
#'   - allocation: Treatment allocation vector (0/1)
#'   - matched_pairs: List of matched pairs with their information
#'   - n_matched_pairs: Number of successfully matched pairs
#'   - n_matched_units: Total number of matched units
#'   - iterations: Number of iterations (always 1 for HCM)
#' 
#' @examples
#' # Generate example data
#' set.seed(2020)
#' n <- 100
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' X <- scale(X)
#' 
#' # Run HCM
#' result <- HCM(X)
#' print(paste("Matched", result$n_matched_pairs, "pairs"))
#' 
#' @export
HCM <- function(X, seed = 2020, method = "single") {
  set.seed(seed)
  n <- nrow(X)
  original_indices <- 1:n
  
  # Check if sample size is even
  if (n %% 2 != 0) {
    warning("Sample size is odd. One unit will be randomly assigned.")
  }
  
  matched_pairs <- list()
  pair_count <- 0
  
  working_data <- X
  working_indices <- original_indices
  
  # Continue until we have fewer than 2 units
  while (nrow(working_data) >= 2) {
    if (nrow(working_data) == 2) {
      # Last two units - create final pair
      pair_count <- pair_count + 1
      treatment_assignment <- sample(c(0, 1), 2, replace = FALSE)
      
      matched_pairs[[pair_count]] <- list(
        indices = working_indices,
        treatments = treatment_assignment,
        distance = as.matrix(dist(working_data))[1, 2]
      )
      break
    }
    
    # Build hierarchical clustering tree
    dist_mat <- dist(working_data, method = "euclidean")
    hc_tree <- hclust(dist_mat, method = method)
    
    # Find the closest pair from the first merge
    first_merge <- hc_tree$merge[1, ]
    
    # Extract indices of the closest pair
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
    
    # Ensure valid indices
    if (idx1 > nrow(working_data)) idx1 <- 1
    if (idx2 > nrow(working_data)) idx2 <- 2
    if (idx1 == idx2) idx2 <- min(idx1 + 1, nrow(working_data))
    
    # Create the pair
    pair_count <- pair_count + 1
    treatment_assignment <- sample(c(0, 1), 2, replace = FALSE)
    pair_distance <- as.matrix(dist_mat)[idx1, idx2]
    
    matched_pairs[[pair_count]] <- list(
      indices = working_indices[c(idx1, idx2)],
      treatments = treatment_assignment,
      distance = pair_distance
    )
    
    # Remove the matched pair from working data
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
  
  # Handle any unmatched units (for odd sample sizes)
  unmatched <- which(is.na(allocation))
  if (length(unmatched) > 0) {
    allocation[unmatched] <- sample(c(0, 1), length(unmatched), replace = TRUE)
  }
  
  return(list(
    allocation = allocation,
    matched_pairs = matched_pairs,
    n_matched_pairs = pair_count,
    n_matched_units = pair_count * 2,
    iterations = 1,
    method = "HCM"
  ))
}


# Helper function to visualize HCM results
#' Visualize HCM Results
#' 
#' Creates diagnostic plots for HCM allocation results
#' 
#' @param X Covariate matrix
#' @param result HCM result object
#' @param max_pairs Maximum number of pairs to display (default: 10)
#' 
#' @export
visualize_HCM <- function(X, result, max_pairs = 10) {
  if (!require(ggplot2)) {
    stop("ggplot2 package required for visualization")
  }
  
  # Create distance distribution plot for matched pairs
  distances <- sapply(result$matched_pairs[1:min(max_pairs, length(result$matched_pairs))], 
                      function(pair) pair$distance)
  
  dist_df <- data.frame(
    pair = 1:length(distances),
    distance = distances
  )
  
  p1 <- ggplot(dist_df, aes(x = pair, y = distance)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    labs(title = "Within-Pair Distances (First 10 Pairs)",
         x = "Pair Number", y = "Euclidean Distance") +
    theme_minimal()
  
  print(p1)
  
  # Treatment allocation summary
  cat("HCM Results Summary:\n")
  cat(paste("Total units:", length(result$allocation), "\n"))
  cat(paste("Matched pairs:", result$n_matched_pairs, "\n"))
  cat(paste("Matched units:", result$n_matched_units, "\n"))
  cat(paste("Treatment group size:", sum(result$allocation), "\n"))
  cat(paste("Control group size:", sum(1 - result$allocation), "\n"))
  cat(paste("Average within-pair distance:", 
            round(mean(sapply(result$matched_pairs, function(p) p$distance)), 4), "\n"))
}