#' Utility Functions for Experimental Design Methods
#' 
#' This file contains helper functions used across different experimental design methods.
#' 
#' @version 1.0.0

# Load required libraries
if (!require("MASS")) install.packages("MASS")
if (!require("stats")) install.packages("stats")

library(MASS)
library(stats)

#' Calculate Mahalanobis Distance for Treatment Allocation
#'
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' @param lambda Ridge parameter (default: 0)
#' @return Mahalanobis distance value
#' @export
maha_dist <- function(X, w, lambda = 0) {
  # Ensure X is centered
  X <- scale(X, center = TRUE, scale = FALSE)
  
  n1 <- sum(w)
  n0 <- sum(1 - w)
  n <- n0 + n1
  
  # Calculate group means
  xbar1 <- t(X) %*% w / n1
  xbar0 <- t(X) %*% (1 - w) / n0
  delta <- xbar1 - xbar0
  
  # Calculate pooled covariance matrix
  cov_mat <- (n / (n1 * n0 * (n - 1))) * t(X) %*% X
  
  # Use generalized inverse for numerical stability
  inv_cov <- MASS::ginv(cov_mat + lambda * diag(nrow(cov_mat)))
  
  # Calculate Mahalanobis distance
  mdist <- as.numeric(t(delta) %*% inv_cov %*% delta)
  return(mdist)
}

#' Calculate Standardized Mean Differences
#'
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' @return Vector of standardized mean differences
#' @export
calc_smd <- function(X, w) {
  treated <- X[w == 1, , drop = FALSE]
  control <- X[w == 0, , drop = FALSE]
  
  mean_t <- colMeans(treated)
  mean_c <- colMeans(control)
  
  var_t <- apply(treated, 2, var)
  var_c <- apply(control, 2, var)
  pooled_sd <- sqrt((var_t + var_c) / 2)
  
  smd_values <- (mean_t - mean_c) / pooled_sd
  return(smd_values)
}

#' Create Greedy Distance-Based Pairs
#'
#' @param X Covariate matrix
#' @return List of paired indices
create_greedy_pairs <- function(X) {
  n <- nrow(X)
  dist_mat <- as.matrix(dist(X))
  
  pairs <- list()
  unmatched <- 1:n
  pair_count <- 0
  
  while (length(unmatched) >= 2) {
    if (length(unmatched) == 2) {
      pair_count <- pair_count + 1
      pairs[[pair_count]] <- unmatched
      break
    }
    
    # Find closest pair among unmatched units
    sub_dist <- dist_mat[unmatched, unmatched]
    diag(sub_dist) <- Inf
    
    min_pos <- which(sub_dist == min(sub_dist), arr.ind = TRUE)[1, ]
    idx1 <- unmatched[min_pos[1]]
    idx2 <- unmatched[min_pos[2]]
    
    pair_count <- pair_count + 1
    pairs[[pair_count]] <- c(idx1, idx2)
    
    unmatched <- unmatched[!unmatched %in% c(idx1, idx2)]
  }
  
  return(pairs)
}

#' Generate Treatment Allocation from Pairs
#'
#' @param pairs List of paired indices
#' @param n Total sample size
#' @return Treatment allocation vector
generate_allocation_from_pairs <- function(pairs, n) {
  allocation <- rep(NA, n)
  
  for (pair in pairs) {
    if (length(pair) == 2) {
      treat_assign <- sample(c(0, 1), 2)
      allocation[pair] <- treat_assign
    }
  }
  
  # Handle any remaining unmatched units
  unmatched <- which(is.na(allocation))
  if (length(unmatched) > 0) {
    allocation[unmatched] <- sample(c(0, 1), length(unmatched), replace = TRUE)
  }
  
  return(allocation)
}
