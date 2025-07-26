#' Hierarchical Clustering Matching with Rerandomization (HCM-ReR)
#' 
#' Combines hierarchical clustering-based pairing with PCA-weighted 
#' rerandomization to achieve superior covariate balance. This method first
#' creates pairs using hierarchical clustering in PCA-weighted space, then
#' applies rerandomization within pairs using a weighted balance criterion.
#' 
#' @param pa Acceptance probability for rerandomization threshold
#' @param X Matrix of covariates (n x p)
#' @param var_explained Proportion of variance explained by PCA components (default: 0.7)
#' @param n_budget Maximum number of rerandomization iterations (default: 1000)
#' @param seed Random seed for reproducibility
#' @param n_sim_threshold Number of simulations for threshold estimation (default: 500)
#' @param clustering_method Hierarchical clustering linkage method (default: "single")
#' 
#' @return List containing:
#'   - allocation: Treatment allocation vector (0/1)
#'   - w: Treatment allocation vector (alias for allocation)
#'   - ii: Number of iterations used
#'   - indicator: TRUE if budget exhausted, FALSE if threshold met
#'   - a: Acceptance threshold used
#'   - k: Number of PCA components used
#'   - matched_pairs: List of matched pairs with their information
#'   - n_matched_pairs: Number of successfully matched pairs
#'   - n_matched_units: Total number of matched units
#'   - iterations: Number of iterations used (alias for ii)
#'   - balance_history: Vector of balance scores across iterations
#'   - method: Method name identifier
#' 
#' @details
#' The HCM-ReR method works in several steps:
#' 1. Perform PCA on covariates and calculate variable importance weights
#' 2. Create hierarchical clustering pairs in PCA-weighted space
#' 3. Estimate acceptance threshold using simulation
#' 4. Perform rerandomization within pairs until threshold is met or budget exhausted
#' 
#' The weighted balance metric gives higher importance to covariates that 
#' contribute more to the first k principal components, allowing for more
#' targeted balance improvement.
#' 
#' @examples
#' # Generate example data with correlation structure
#' set.seed(2020)
#' n <- 200
#' p <- 20
#' Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))  # AR(1) correlation
#' X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
#' X <- scale(X)
#' 
#' # Run HCM-ReR with default settings
#' result <- HCM_ReR(pa = 0.05, X = X)
#' 
#' # Check results
#' print(paste("Converged in", result$iterations, "iterations"))
#' print(paste("Final balance score:", round(tail(result$balance_history, 1), 4)))
#' 
#' # Run with custom settings
#' result_custom <- HCM_ReR(pa = 0.1, X = X, var_explained = 0.8, n_budget = 2000)
#' 
#' @export
HCM_ReR <- function(pa, X, var_explained = 0.7, n_budget = 1000, 
                    seed = 2020, n_sim_threshold = 500, 
                    clustering_method = "single") {
  
  # Input validation
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (pa <= 0 || pa >= 1) {
    stop("pa must be between 0 and 1")
  }
  if (n_budget <= 0) {
    stop("n_budget must be positive")
  }
  
  set.seed(seed)
  n <- nrow(X)
  d <- ncol(X)
  
  if (n %% 2 != 0) {
    stop("Sample size must be even for pair-based methods")
  }
  
  # Step 1: Calculate PCA weights for variable importance
  X_svd <- svd(X)
  
  if (var_explained == 'Kaiser') {
    # Kaiser criterion: keep components with eigenvalues > mean eigenvalue
    k <- sum(X_svd$d^2 > mean(X_svd$d^2))
  } else {
    # Cumulative variance explained criterion
    cumsum_vars <- cumsum(X_svd$d^2)
    cumsum_vars <- cumsum_vars / cumsum_vars[length(cumsum_vars)]
    k <- sum(cumsum_vars < var_explained) + 1
  }
  
  # Ensure we have at least 1 component
  k <- max(1, min(k, d))
  
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
    
    # Standardize covariates for fair comparison
    X_std <- scale(X)
    mean_diff <- colMeans(X_std[W == 1, , drop = FALSE]) - 
      colMeans(X_std[W == 0, , drop = FALSE])
    
    # Weighted sum of squared standardized mean differences
    weighted_balance <- sum(weights * mean_diff^2)
    return(weighted_balance)
  }
  
  # Step 3: Estimate acceptance threshold via simulation
  sim_distances <- replicate(n_sim_threshold, {
    W_sim <- sample(c(rep(1, n/2), rep(0, n/2)))
    calc_weighted_balance(X, W_sim, var_weights)
  })
  threshold <- quantile(sim_distances, pa, na.rm = TRUE)
  
  # Step 4: Create hierarchical clustering pairs in weighted space
  X_std <- scale(X)
  weighted_X_std <- sweep(X_std, 2, sqrt(var_weights), "*")
  
  # Hierarchical clustering in weighted space
  dist_mat <- dist(weighted_X_std, method = "euclidean")
  hc_tree <- hclust(dist_mat, method = clustering_method)
  
  # Extract pairs from hierarchical clustering tree
  pairs <- list()
  used_indices <- logical(n)
  pair_count <- 0
  
  # Process merges to create pairs
  for (i in 1:nrow(hc_tree$merge)) {
    merge_row <- hc_tree$merge[i, ]
    
    # Get actual indices (negative values indicate original observations)
    idx1 <- if (merge_row[1] < 0) abs(merge_row[1]) else NULL
    idx2 <- if (merge_row[2] < 0) abs(merge_row[2]) else NULL
    
    # If both are original points and neither is used, create a pair
    if (!is.null(idx1) && !is.null(idx2) && 
        !used_indices[idx1] && !used_indices[idx2]) {
      pair_count <- pair_count + 1
      pairs[[pair_count]] <- c(idx1, idx2)
      used_indices[c(idx1, idx2)] <- TRUE
      
      # Stop when we have enough pairs
      if (pair_count >= n/2) break
    }
  }
  
  # Handle any remaining unpaired units with greedy pairing
  remaining <- which(!used_indices)
  if (length(remaining) >= 2) {
    remaining_pairs <- create_greedy_pairs(weighted_X_std[remaining, , drop = FALSE])
    for (i in seq_along(remaining_pairs)) {
      pair_count <- pair_count + 1
      pairs[[pair_count]] <- remaining[remaining_pairs[[i]]]
    }
  }
  
  # Step 5: Rerandomization within pairs using weighted balance criterion
  best_w <- NULL
  best_balance <- Inf
  balance_history <- numeric()
  
  for (iter in 1:n_budget) {
    # Generate allocation within current pairs
    w <- generate_allocation_from_pairs(pairs, n)
    
    # Skip if allocation is invalid
    if (any(is.na(w)) || sum(w) == 0 || sum(w) == length(w)) {
      balance_score <- Inf
    } else {
      # Evaluate using weighted balance metric
      balance_score <- calc_weighted_balance(X, w, var_weights)
    }
    
    balance_history <- c(balance_history, balance_score)
    
    # Update best allocation found so far
    if (balance_score < best_balance) {
      best_w <- w
      best_balance <- balance_score
    }
    
    # Check if threshold is met
    if (balance_score <= threshold) {
      matched_pairs_info <- convert_pairs_to_matched_format(pairs, w, X)
      return(list(
        allocation = w, 
        w = w, 
        ii = iter, 
        indicator = FALSE, 
        a = threshold, 
        k = k,
        matched_pairs = matched_pairs_info$matched_pairs,
        n_matched_pairs = matched_pairs_info$n_matched_pairs,
        n_matched_units = matched_pairs_info$n_matched_units,
        iterations = iter, 
        balance_history = balance_history,
        method = "HCM-ReR"
      ))
    }
  }
  
  # Return best allocation if threshold not met within budget
  if (is.null(best_w)) {
    best_w <- sample(c(rep(1, n/2), rep(0, n/2)))
  }
  
  matched_pairs_info <- convert_pairs_to_matched_format(pairs, best_w, X)
  return(list(
    allocation = best_w, 
    w = best_w, 
    ii = n_budget, 
    indicator = TRUE, 
    a = threshold, 
    k = k,
    matched_pairs = matched_pairs_info$matched_pairs,
    n_matched_pairs = matched_pairs_info$n_matched_pairs,
    n_matched_units = matched_pairs_info$n_matched_units,
    iterations = n_budget, 
    balance_history = balance_history,
    method = "HCM-ReR"
  ))
}


# Helper functions for HCM-ReR
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

generate_allocation_from_pairs <- function(pairs, n) {
  allocation <- rep(NA, n)
  
  for (pair in pairs) {
    if (length(pair) == 2) {
      treat_assign <- sample(c(0, 1), 2)
      allocation[pair] <- treat_assign
    }
  }
  
  return(allocation)
}

convert_pairs_to_matched_format <- function(pairs, allocation, X) {
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
    matched_pairs = matched_pairs,
    n_matched_pairs = length(matched_pairs),
    n_matched_units = length(matched_pairs) * 2
  ))
}


#' Visualize HCM-ReR Convergence
#' 
#' Creates diagnostic plots for HCM-ReR rerandomization process
#' 
#' @param result HCM-ReR result object
#' @param show_threshold Whether to show acceptance threshold line
#' 
#' @export
visualize_HCM_ReR_convergence <- function(result, show_threshold = TRUE) {
  if (!require(ggplot2)) {
    stop("ggplot2 package required for visualization")
  }
  
  balance_df <- data.frame(
    iteration = 1:length(result$balance_history),
    balance_score = result$balance_history
  )
  
  p <- ggplot(balance_df, aes(x = iteration, y = balance_score)) +
    geom_line(color = "blue", alpha = 0.7) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(title = "HCM-ReR Convergence",
         subtitle = paste("Converged in", result$iterations, "iterations"),
         x = "Iteration", y = "Weighted Balance Score") +
    theme_minimal()
  
  if (show_threshold && !is.null(result$a)) {
    p <- p + geom_hline(yintercept = result$a, color = "red", linetype = "dashed",
                        alpha = 0.7) +
      annotate("text", x = length(result$balance_history) * 0.8, y = result$a,
               label = paste("Threshold =", round(result$a, 4)), 
               vjust = -0.5, color = "red")
  }
  
  print(p)
  
  # Print summary statistics
  cat("HCM-ReR Convergence Summary:\n")
  cat(paste("Total iterations:", result$iterations, "\n"))
  cat(paste("Budget exhausted:", result$indicator, "\n"))
  cat(paste("Final balance score:", round(tail(result$balance_history, 1), 6), "\n"))
  cat(paste("Acceptance threshold:", round(result$a, 6), "\n"))
  cat(paste("PCA components used:", result$k, "\n"))
  cat(paste("Improvement over iteration 1:", 
            round((result$balance_history[1] - tail(result$balance_history, 1)) / 
                    result$balance_history[1] * 100, 2), "%\n"))
}