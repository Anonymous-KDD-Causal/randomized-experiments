# ============================================================================
# RANDOMIZATION METHODS COMPARISON ON REAL-WORLD DATA
# ============================================================================

library(glmnet)
library(fastDummies)
library(latex2exp)
library(ggplot2)
library(MASS)
library(lattice)
library(reshape2)
library(scales)

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

# Mahalanobis distance calculation
maha_dist = function(X,w,lambda=0){
  # note that the input data matrix should be centered
  library(MASS)
  n1 = sum(w)
  n0 = sum(1-w)
  n = n0+n1
  xbar1 = t(X)%*%w/n1
  xbar0 = t(X)%*%(1-w)/n0
  delta = xbar1-xbar0
  cov_mat = (n/(n1*n0*(n-1)))*t(X)%*%X
  # use ginv other than solve to prevent some numerical issues
  mdist = as.numeric(t(delta)%*%ginv(cov_mat+lambda*diag(nrow(cov_mat)))%*%delta)
  return(mdist)
}

# Helper functions
mean_diff <- function(X, w) {
  n1 <- sum(w)
  n0 <- sum(1-w)
  xbar1 <- t(X) %*% w / n1
  xbar0 <- t(X) %*% (1-w) / n0
  delta <- xbar1 - xbar0
  return(delta)
}

tau_est <- function(y, w) {
  mean(y[w==1]) - mean(y[w==0])
}

# Standardized mean difference calculation
smd <- function(X, w) {
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

# Function to create implicit matched pairs for balanced allocations
create_implicit_pairs <- function(X, allocation, method = "random") {
  treated_indices <- which(allocation == 1)
  control_indices <- which(allocation == 0)
  
  n_treated <- length(treated_indices)
  n_control <- length(control_indices)
  n_pairs <- min(n_treated, n_control)
  
  matched_pairs <- list()
  
  if (method == "random") {
    # Random pairing
    if (n_pairs > 0) {
      treated_sample <- sample(treated_indices, n_pairs)
      control_sample <- sample(control_indices, n_pairs)
    } else {
      treated_sample <- control_sample <- integer(0)
    }
  } else if (method == "distance") {
    # Distance-based pairing (greedy nearest neighbor)
    if (n_pairs > 0) {
      treated_sample <- treated_indices[1:n_pairs]
      control_sample <- control_indices[1:n_pairs]
      
      # Calculate distances and do greedy matching
      dist_matrix <- as.matrix(dist(X[c(treated_sample, control_sample), ]))
      n_t <- length(treated_sample)
      
      used_control <- rep(FALSE, n_pairs)
      final_control <- numeric(n_pairs)
      
      for (i in 1:n_pairs) {
        available_control <- which(!used_control)
        if (length(available_control) == 0) break
        
        # Find closest available control for this treated unit
        distances_to_controls <- dist_matrix[i, (n_t + 1):(n_t + n_pairs)]
        distances_to_controls[used_control] <- Inf
        
        best_control_idx <- which.min(distances_to_controls)
        final_control[i] <- best_control_idx
        used_control[best_control_idx] <- TRUE
      }
      
      control_sample <- control_indices[final_control]
    } else {
      treated_sample <- control_sample <- integer(0)
    }
  }
  
  # Create matched pairs
  if (n_pairs > 0) {
    for (i in 1:n_pairs) {
      pair_distance <- dist(X[c(treated_sample[i], control_sample[i]), ])
      
      matched_pairs[[i]] <- list(
        indices = c(treated_sample[i], control_sample[i]),
        treatments = c(1, 0),
        distance = as.numeric(pair_distance)
      )
    }
  }
  
  return(list(
    matched_pairs = matched_pairs,
    n_matched_pairs = n_pairs,
    n_matched_units = n_pairs * 2
  ))
}

# ============================================================================
# HELPER FUNCTIONS FOR PAIR-BASED METHODS
# ============================================================================

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

# ============================================================================
# HIERARCHICAL CLUSTERING MATCHING (HCM)
# ============================================================================

hcm_randomized <- function(X, linkage = "single", seed = 2020) {
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
    iterations = 1
  ))
}

# ============================================================================
# HCM-ReR METHOD
# ============================================================================

HCM_ReR <- function(X, pa = 0.2, n_budget = 1000, linkage = "single", seed = 2020) {
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
  
  library(parallel)
  library(MASS)
  
  # Windows-compatible parallel processing
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  
  # Export all necessary functions and objects
  clusterExport(cl, c("generate_allocation_from_pairs", "maha_dist", "pairs", "X", "n",
                      "create_greedy_pairs"))  # Add any other helper functions used
  
  # Load required libraries on each worker
  clusterEvalQ(cl, {
    library(MASS)
  })
  
  sim_distances <- parSapply(cl, 1:400, function(i) {
    allocation_sim <- generate_allocation_from_pairs(pairs, n)
    maha_dist(X, allocation_sim)
  })
  
  stopCluster(cl)
  threshold <- quantile(sim_distances, pa, na.rm = TRUE)
  
  # Step 2: Rerandomization within pairs
  best_allocation <- NULL
  best_balance <- Inf
  iterations_used <- n_budget
  
  for (iter in 1:n_budget) {
    # Generate allocation from pairs
    print(iter)
    allocation <- generate_allocation_from_pairs(pairs, n)
    
    # Calculate balance
    current_balance <- maha_dist(X, allocation)
    
    if (current_balance < best_balance) {
      best_allocation <- allocation
      best_balance <- current_balance
    }
    
    # Check if threshold is met
    if (current_balance <= threshold) {
      iterations_used <- iter
      best_allocation <- allocation
      break
    }
  }
  
  # Convert pairs to matched pair format using final allocation
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
    n_matched_pairs = pair_count,
    n_matched_units = pair_count * 2,
    iterations = iterations_used
  ))
}
# ============================================================================
# RERANDOMIZATION METHODS
# ============================================================================

# Standard Rerandomization
ReR <- function(pa, X, n_budget=1000, seed=2020) {
  set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  a <- qchisq(pa, d)
  
  n0 <- ceiling(n/2)
  n1 <- n - n0
  w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
  mdist <- ifelse(n < d, n-1, maha_dist(X, w))
  
  ii <- 1
  best_w <- w
  best_mdist <- mdist
  while (mdist > a) {
    w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
    mdist <- maha_dist(X, w)
    if(best_mdist > mdist){
      best_w <- w
      best_mdist <- mdist
    }
    ii <- ii + 1
    if(ii >= n_budget){
      break
    }
  }
  
  indicator <- ii >= n_budget
  if(indicator == TRUE){
    return_w <- best_w
  } else {
    return_w <- w
  }
  
  pair_info <- create_implicit_pairs(X, return_w, method = "distance")
  
  return(list(
    allocation = return_w,
    w = return_w, 
    ii = ii, 
    indicator = indicator, 
    a = a,
    matched_pairs = pair_info$matched_pairs,
    n_matched_pairs = pair_info$n_matched_pairs,
    n_matched_units = pair_info$n_matched_units,
    iterations = ii
  ))
}

# PCA Rerandomization
PCA_ReR <- function(pa, X, var_ratio=0.7, n_budget=1000, seed=2020) {
  set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  
  X_svd <- svd(X)
  if (var_ratio == 'Kaiser'){
    k <- sum(X_svd$d^2 > mean(X_svd$d^2))
  } else {
    cumsum_vars <- cumsum(X_svd$d^2)
    cumsum_vars <- cumsum_vars / cumsum_vars[length(cumsum_vars)]
    k <- sum(cumsum_vars < var_ratio) + 1
  }
  
  if(k == 1){
    Zk <- as.matrix(X_svd$u[,1] * X_svd$d[1])
  } else {
    Zk <- X_svd$u[,1:k] %*% diag(X_svd$d[1:k])
  }
  a_pca <- qchisq(pa, k)
  
  n0 <- ceiling(n/2)
  n1 <- n - n0
  w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
  mdist <- maha_dist(Zk, w)
  
  ii <- 1
  best_w <- w
  best_mdist <- mdist
  while (mdist > a_pca) {
    w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
    mdist <- maha_dist(Zk, w)
    if(best_mdist > mdist){
      best_w <- w
      best_mdist <- mdist
    }
    ii <- ii + 1
    if(ii >= n_budget){
      break
    }
  }
  
  indicator <- ii >= n_budget
  if(indicator == TRUE){
    return_w <- best_w
  } else {
    return_w <- w
  }
  
  pair_info <- create_implicit_pairs(X, return_w, method = "distance")
  
  return(list(
    allocation = return_w,
    k = k, 
    w = return_w, 
    indicator = indicator, 
    ii = ii, 
    a = a_pca,
    matched_pairs = pair_info$matched_pairs,
    n_matched_pairs = pair_info$n_matched_pairs,
    n_matched_units = pair_info$n_matched_units,
    iterations = ii
  ))
}

# Ridge Rerandomization helper functions
find_q <- function(lambda, lambdas, pa=0.05, ub=NULL) {
  d <- length(lambdas)
  f <- function(t, q){
    tmp1 <- sin(0.5*(-t*q + sum(atan(lambdas/(lambda+lambdas+1e-8)*t))))
    tmp2 <- t*prod(1+(lambdas/(lambda+lambdas+1e-8))^2*t^2)^(1/4)
    tmp1/tmp2
  }
  vf <- Vectorize(f, vectorize.args = 't')
  Fq <- function(q){
    tryCatch({
      0.5-(1/pi)*integrate(vf, lower = 0, upper = Inf, q=q, subdivisions=5000)$value
    }, error=function(e){
      xi <- 1e-4
      U <- (xi*pi*(d/2)*prod(lambdas/(lambda+lambdas+1e-8)))^(2/d)
      0.5-(1/pi)*integrate(vf, lower = 0, upper = U, q=q, subdivisions=5000)$value
    })
  }
  fn <- function(q, pa){
    (Fq(q)-pa)^2
  }
  
  seq_ub <- ifelse(is.null(ub), d, ub)
  if(seq_ub > 0.2){
    init_seq <- seq(0.1, seq_ub, 0.1)
  } else {
    init_seq <- seq(0, seq_ub, length=100)
  }
  init_idx <- which.min((sapply(init_seq, Fq)-pa)^2)
  result <- nlminb(start=init_seq[init_idx], objective=fn, lower=0, pa=pa)$par
  result
}

find_lambda <- function(sigmat, lambdas, V, pa=0.05, delta=1e-2, eps=1e-4, n_search=10, nn=1000) {
  d <- nrow(V)
  lambda <- 0 
  Lambda <- c()
  dk_mat <- c()
  a_seq <- c()
  Z <- matrix(rnorm(nn*d), nn, d)
  
  a <- find_q(lambda, lambdas, pa)
  a_plus <- find_q(lambda+delta, lambdas, pa, a)
  
  i <- 0
  while (abs((lambda+delta)*a_plus-lambda*a) > eps) {
    lambda <- lambda + delta
    a <- a_plus
    a_plus <- find_q(lambda+delta, lambdas, pa, a)
    
    M_lambda <- Z^2 %*% (lambdas/(lambda+lambdas))
    indicator <- M_lambda <= a
    dk <- t(indicator) %*% Z^2 / (sum(indicator)+eps)
    vk <- diag(V %*% diag(lambdas*as.vector(dk)) %*% t(V)) / (diag(sigmat)+1e-8)
    if(mean(vk) < pchisq(qchisq(pa, df=d), d+2)/pa){
      Lambda <- c(Lambda, lambda)
      dk_mat <- rbind(dk_mat, dk)
      a_seq <- c(a_seq, a)
    }
    
    i <- i + 1
    if(i > n_search){
      break
    }
  }
  
  if(length(Lambda) == 0){
    lambda <- 0
    a <- find_q(0, lambdas, pa)
  } else {
    ck <- lambdas^2 / sum(lambdas^2)
    objs <- dk_mat^2 %*% ck - (dk_mat %*% ck)^2
    lambda <- Lambda[which.min(objs)]
    a <- a_seq[which.min(objs)]
  }
  
  return(list(lambda=lambda, a=a))
}

# Ridge Rerandomization
RidgeReR <- function(pa, X, n_budget=1000, seed=2020) {
  set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  
  sigmat <- cov(X) * 4/n
  eig_decomp <- eigen(sigmat)
  lambdas <- eig_decomp$values
  V <- eig_decomp$vectors
  
  params <- find_lambda(sigmat, lambdas, V, pa)
  a <- params$a
  lambda <- params$lambda
  
  n0 <- ceiling(n/2)
  n1 <- n - n0
  w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
  mdist <- maha_dist(X, w, lambda)
  
  ii <- 1
  best_w <- w
  best_mdist <- mdist
  while (mdist > a) {
    w <- sample(rep(c(0,1), c(n0,n1)), n, FALSE)
    mdist <- maha_dist(X, w, lambda)
    if(best_mdist > mdist){
      best_w <- w
      best_mdist <- mdist
    }
    ii <- ii + 1
    if(ii >= n_budget){
      break
    }
  }
  
  indicator <- ii >= n_budget
  if(indicator == TRUE){
    return_w <- best_w
  } else {
    return_w <- w
  }
  
  pair_info <- create_implicit_pairs(X, return_w, method = "distance")
  
  return(list(
    allocation = return_w,
    w = return_w, 
    indicator = indicator, 
    ii = ii, 
    a = a,
    matched_pairs = pair_info$matched_pairs,
    n_matched_pairs = pair_info$n_matched_pairs,
    n_matched_units = pair_info$n_matched_units,
    iterations = ii
  ))
}

# Simple Random Allocation
simple_randomization <- function(X, seed = 2020) {
  set.seed(seed)
  n <- nrow(X)
  allocation <- sample(c(rep(0, n/2), rep(1, n/2)), n, replace = FALSE)
  
  pair_info <- create_implicit_pairs(X, allocation, method = "random")
  
  return(list(
    allocation = allocation,
    matched_pairs = pair_info$matched_pairs,
    n_matched_pairs = pair_info$n_matched_pairs,
    n_matched_units = pair_info$n_matched_units,
    iterations = 1
  ))
}

# Corrected PCA-Weighted Distance Rerandomization (PWD_ReR)
PWD_ReR <- function(pa, X, var_explained = 0.7, n_budget = 1000, 
                    seed = 2020, n_sim_threshold=500) {
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  n0 <- ceiling(n/2)
  n1 <- n - n0
  
  # Calculate PCA weights (same as before)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_svd <- svd(X_centered)
  
  if (var_explained == 'Kaiser') {
    k <- sum(X_svd$d^2 > mean(X_svd$d^2))
  } else {
    cumsum_vars <- cumsum(X_svd$d^2)
    cumsum_vars <- cumsum_vars / cumsum_vars[length(cumsum_vars)]
    k <- sum(cumsum_vars < var_explained) + 1
  }
  
  rotation <- X_svd$v[, 1:k, drop = FALSE]
  
  if (k == 1) {
    var_weights <- rotation^2
  } else {
    var_weights <- rowSums(rotation^2)
  }
  var_weights <- var_weights / sum(var_weights)
  # Simple weighted sum of squared standardized mean differences
  calc_weighted_balance <- function(X, W, weights) {
    n1 <- sum(W)
    n0 <- sum(1-W)
    if (n1 == 0 || n0 == 0) return(Inf)
    
    # Standardize X for fair comparison across variables
    X_std <- scale(X)
    mean_diff <- colMeans(X_std[W == 1, , drop = FALSE]) - colMeans(X_std[W == 0, , drop = FALSE])
    
    # Weighted sum of squared differences
    weighted_balance <- sum(weights * mean_diff^2)
    return(weighted_balance)
  }
  
  # Simulate correct threshold for weighted distance
  sim_distances <- replicate(n_sim_threshold, {
    W_sim <- sample(c(rep(1, n1), rep(0, n0)))
    calc_weighted_balance(X, W_sim, var_weights)
  })
  
  threshold <- quantile(sim_distances, pa, na.rm = TRUE)
  
  best_balance <- Inf
  best_W <- NULL
  
  for(iter in 1:n_budget) {
    W <- sample(c(rep(1, n1), rep(0, n0)))
    
    balance_score <- calc_weighted_balance(X, W, var_weights)
    
    if(balance_score < best_balance) {
      best_balance <- balance_score
      best_W <- W
    }
    
    if(balance_score <= threshold) {
      # Accept this allocation
      pair_info <- create_implicit_pairs(X, W, method = "distance")
      return(list(
        allocation = W,
        w = W,
        indicator = FALSE,
        ii = iter,
        k = k,
        a = threshold,
        matched_pairs = pair_info$matched_pairs,
        n_matched_pairs = pair_info$n_matched_pairs,
        n_matched_units = pair_info$n_matched_units,
        iterations = iter
      ))
    }
  }
  
  # Return best if no allocation accepted
  pair_info <- create_implicit_pairs(X, best_W, method = "distance")
  return(list(
    allocation = best_W,
    w = best_W,
    indicator = TRUE,
    ii = n_budget,
    k = k,
    a = threshold,
    matched_pairs = pair_info$matched_pairs,
    n_matched_pairs = pair_info$n_matched_pairs,
    n_matched_units = pair_info$n_matched_units,
    iterations = n_budget
  ))
}
# ============================================================================
# REAL IHDP DATA LOADING AND PREPROCESSING
# Based on the provided CEVAE repository approach
# ============================================================================

# load functions
src_dir = 'R/data/'
files = list.files(src_dir)
for(file in files){
  source(paste(src_dir,file,sep=''))
}

# load the data from Hill (2011)
load('R/data/example.data')
load('R/data/sim.data')

# combine imp1 with the treatment effect
imp1 = data.frame(iqsb.36=ihdp$iqsb.36,imp1)
covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
p=length(c(covs.cont.n,covs.cat.n))
imp1 = imp1[,c('iqsb.36','treat',covs.cont.n,covs.cat.n)]

for(cov_name in covs.cat.n){
  cat(cov_name,unique(imp1[,cov_name]),'\n')
}
imp1[,'first'] = ifelse(imp1[,'first']==2,1,0)
imp1=na.omit(imp1) # removing the missing values
ihdp = imp1

# Data Preprocessing ------------------------------------------------------

# this data set can be found on the repo
# https://github.com/AMLab-Amsterdam/CEVAE

colnames(ihdp) = c('y','treatment',
                   paste('x',1:(ncol(ihdp)-2),sep=''))
# standardize all covariates
ihdp[,1:27] = scale(ihdp[,1:27])

# Simulation --------------------------------------------------------------
# fit a linear model to simulate the real process
# all variables
set.seed(2020)
vars = paste(paste('x',1:(ncol(ihdp)-2),sep=''),collapse = '+')
# formula
fml = as.formula(paste('~ (',vars,') ^2 - 1'))
design_mat = data.frame(treatment=ihdp$treatment,
                        model.matrix(fml,data=ihdp))
design_mat = scale(design_mat)
cv_output = cv.glmnet(x=as.matrix(design_mat),y=ihdp$y,
                      family = 'gaussian')
best_lam <- cv_output$lambda.min
lm_fit = glmnet(x=as.matrix(design_mat),y=ihdp$y,
                lambda=best_lam,
                family = 'gaussian')
tau = coef(lm_fit)[2]

# a function to generate the response
ihdp_y_gen = function(X,w){
  XX = cbind(w,X)
  yhat = predict(lm_fit,XX)
  yhat
}

# We balance all samples
n = nrow(ihdp)
X = as.matrix(model.matrix(fml,data=ihdp))[1:n,]
X = scale(X)

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

reps <- 1000
n_budget <- 1000
pa <- 0.05
var_ratio <- 0.7

# Create save directory
save_nm <- 'selected_methods_study'
save_dir_folder <- file.path('R/data/', save_nm)
if (!file.exists(save_dir_folder)) {
  dir.create(save_dir_folder, recursive = TRUE)
}

# ============================================================================
# RUN SELECTED METHODS ONLY
# ============================================================================

cat("Starting Selected Methods Simulation Study...\n")
cat("Sample size:", n, "\n")
cat("Number of covariates:", ncol(X), "\n")
cat("Number of replications:", reps, "\n")
cat("Budget per replication:", n_budget, "\n\n")

# Method 1: Simple Randomization (Baseline)
cat("Running Simple Randomization...\n")
save_dir_simple <- file.path(save_dir_folder, 'simple.rdata')
if (!file.exists(save_dir_simple)) {
  set.seed(2020)
  ii_vec_simple <- c()
  W_simple <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- simple_randomization(X, seed = 2020 + i)
    ii_vec_simple <- c(ii_vec_simple, result$iterations)
    W_simple <- cbind(W_simple, as.numeric(result$allocation))
    if (i %% 100 == 0) cat("Completed", i, "replications\n")
  }
  end_time <- Sys.time()
  time_simple <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_simple, time_simple, ii_vec_simple, file = save_dir_simple)
} else {
  load(file = save_dir_simple)
}

# Method 2: Standard Rerandomization
cat("Running Standard Rerandomization...\n")
save_dir_rer <- file.path(save_dir_folder, 'rer.rdata')
if (!file.exists(save_dir_rer)) {
  set.seed(2020)
  ii_vec_rer <- c()
  W_rer <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- ReR(pa = pa, X = X, n_budget = n_budget, seed = 2020 + i)
    ii_vec_rer <- c(ii_vec_rer, result$ii)
    W_rer <- cbind(W_rer, as.numeric(result$w))
    if (i %% 100 == 0) cat("Completed", i, "replications\n")
  }
  end_time <- Sys.time()
  time_rer <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_rer, time_rer, ii_vec_rer, file = save_dir_rer)
} else {
  load(file = save_dir_rer)
}

# Method 3: Ridge Rerandomization
cat("Running Ridge Rerandomization...\n")
save_dir_ridgerer <- file.path(save_dir_folder, 'ridgerer.rdata')
if (!file.exists(save_dir_ridgerer)) {
  set.seed(2020)
  ii_vec_ridgerer <- c()
  W_ridgerer <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- RidgeReR(pa = pa, X = X, n_budget = n_budget, seed = 2020 + i)
    ii_vec_ridgerer <- c(ii_vec_ridgerer, result$ii)
    W_ridgerer <- cbind(W_ridgerer, as.numeric(result$w))
    if (i %% 100 == 0) cat("Completed", i, "replications\n")
  }
  end_time <- Sys.time()
  time_ridgerer <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_ridgerer, time_ridgerer, ii_vec_ridgerer, file = save_dir_ridgerer)
} else {
  load(file = save_dir_ridgerer)
}

# Method 4: PCA Rerandomization
cat("Running PCA Rerandomization...\n")
save_dir_pcarer <- file.path(save_dir_folder, 'pcarer.rdata')
if (!file.exists(save_dir_pcarer)) {
  set.seed(2020)
  ii_vec_pcarer <- c()
  W_pcarer <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- PCA_ReR(pa = pa, X = X, var_ratio = var_ratio, 
                      n_budget = n_budget, seed = 2020 + i)
    ii_vec_pcarer <- c(ii_vec_pcarer, result$ii)
    W_pcarer <- cbind(W_pcarer, as.numeric(result$w))
    if (i %% 100 == 0) cat("Completed", i, "replications\n")
  }
  end_time <- Sys.time()
  time_pcarer <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_pcarer, time_pcarer, ii_vec_pcarer, file = save_dir_pcarer)
} else {
  load(file = save_dir_pcarer)
}

# Method 5: PWD Rerandomization
cat("Running PWD Rerandomization...\n")
save_dir_pwdrer <- file.path(save_dir_folder, 'pwdrer.rdata')
if (!file.exists(save_dir_pwdrer)) {
  set.seed(2020)
  ii_vec_pwdrer <- c()
  W_pwdrer <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- PWD_ReR(pa = pa, X = X, var_explained = var_ratio,
                      n_budget = n_budget, seed = 2020 + i)
    ii_vec_pwdrer <- c(ii_vec_pwdrer, result$ii)
    W_pwdrer <- cbind(W_pwdrer, as.numeric(result$w))
    if (i %% 100 == 0) cat("Completed", i, "replications\n")
  }
  end_time <- Sys.time()
  time_pwdrer <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_pwdrer, time_pwdrer, ii_vec_pwdrer, file = save_dir_pwdrer)
} else {
  load(file = save_dir_pwdrer)
}

# Method 6: Hierarchical Clustering Matching
cat("Running Hierarchical Clustering Matching...\n")
save_dir_hcm <- file.path(save_dir_folder, 'hcm.rdata')
if (!file.exists(save_dir_hcm)) {
  set.seed(2020)
  ii_vec_hcm <- c()
  W_hcm <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- hcm_randomized(X, seed = 2020 + i)
    ii_vec_hcm <- c(ii_vec_hcm, result$iterations)
    W_hcm <- cbind(W_hcm, as.numeric(result$allocation))
    if (i %% 10 == 0) cat("Completed", i, "replications\n")
  }
  end_time <- Sys.time()
  time_hcm <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_hcm, time_hcm, ii_vec_hcm, file = save_dir_hcm)
} else {
  load(file = save_dir_hcm)
}

# Method 7: HCM-ReR
cat("Running HCM-ReR...\n")
save_dir_hcm_rer <- file.path(save_dir_folder, 'hcm_rer.rdata')
if (!file.exists(save_dir_hcm_rer)) {
  set.seed(2020)
  ii_vec_hcm_rer <- c()
  W_hcm_rer <- c()
  
  start_time <- Sys.time()
  for (i in 1:reps) {
    result <- HCM_ReR(pa = pa, X = X, n_budget = n_budget, seed = 2020 + i)
    
    # Fix: Use correct element names
    ii_vec_hcm_rer <- c(ii_vec_hcm_rer, result$iterations)  # Changed from result$ii
    W_hcm_rer <- cbind(W_hcm_rer, as.numeric(result$allocation))  # Changed from result$w
    
    if (i %% 1 == 0) cat("Completed", i, "replications\n")  # Changed to every 100 for less output
  }
  end_time <- Sys.time()
  time_hcm_rer <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  save(W_hcm_rer, time_hcm_rer, ii_vec_hcm_rer, file = save_dir_hcm_rer)
} else {
  load(file = save_dir_hcm_rer)
}

# ============================================================================
# ANALYSIS AND RESULTS
# ============================================================================

cat("\n=== ANALYSIS PHASE ===\n")

# Calculate covariate differences for selected methods
cov_diff_simple <- apply(W_simple, MARGIN = 2, FUN = mean_diff, X = X)
cov_diff_rer <- apply(W_rer, MARGIN = 2, FUN = mean_diff, X = X)
cov_diff_ridgerer <- apply(W_ridgerer, MARGIN = 2, FUN = mean_diff, X = X)
cov_diff_pcarer <- apply(W_pcarer, MARGIN = 2, FUN = mean_diff, X = X)
cov_diff_pwdrer <- apply(W_pwdrer, MARGIN = 2, FUN = mean_diff, X = X)
cov_diff_hcm <- apply(W_hcm, MARGIN = 2, FUN = mean_diff, X = X)
cov_diff_hcm_rer <- apply(W_hcm_rer, MARGIN = 2, FUN = mean_diff, X = X)

# Calculate treatment effects
tau_true <- tau  # True treatment effect from our simulation

tau_simple <- apply(W_simple, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

tau_rer <- apply(W_rer, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

tau_ridgerer <- apply(W_ridgerer, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

tau_pcarer <- apply(W_pcarer, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

tau_pwdrer <- apply(W_pwdrer, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

tau_hcm <- apply(W_hcm, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

tau_hcm_rer <- apply(W_hcm_rer, MARGIN = 2, FUN = function(w) {
  y <- ihdp_y_gen(X, w)
  tau_est(y, w)
})

# Calculate Mahalanobis distances
maha_dist_simple <- apply(W_simple, MARGIN = 2, FUN = function(w) maha_dist(X, w))
maha_dist_rer <- apply(W_rer, MARGIN = 2, FUN = function(w) maha_dist(X, w))
maha_dist_ridgerer <- apply(W_ridgerer, MARGIN = 2, FUN = function(w) maha_dist(X, w))
maha_dist_pcarer <- apply(W_pcarer, MARGIN = 2, FUN = function(w) maha_dist(X, w))
maha_dist_pwdrer <- apply(W_pwdrer, MARGIN = 2, FUN = function(w) maha_dist(X, w))
maha_dist_hcm <- apply(W_hcm, MARGIN = 2, FUN = function(w) maha_dist(X, w))
maha_dist_hcm_rer <- apply(W_hcm_rer, MARGIN = 2, FUN = function(w) maha_dist(X, w))

# Calculate bias and RMSE
bias_simple <- mean(tau_simple - tau_true)
bias_rer <- mean(tau_rer - tau_true)
bias_ridgerer <- mean(tau_ridgerer - tau_true)
bias_pcarer <- mean(tau_pcarer - tau_true)
bias_pwdrer <- mean(tau_pwdrer - tau_true)
bias_hcm <- mean(tau_hcm - tau_true)
bias_hcm_rer <- mean(tau_hcm_rer - tau_true)

rmse_simple <- sqrt(mean((tau_simple - tau_true)^2))
rmse_rer <- sqrt(mean((tau_rer - tau_true)^2))
rmse_ridgerer <- sqrt(mean((tau_ridgerer - tau_true)^2))
rmse_pcarer <- sqrt(mean((tau_pcarer - tau_true)^2))
rmse_pwdrer <- sqrt(mean((tau_pwdrer - tau_true)^2))
rmse_hcm <- sqrt(mean((tau_hcm - tau_true)^2))
rmse_hcm_rer <- sqrt(mean((tau_hcm_rer - tau_true)^2))

# Calculate mean SMD - CORRECTED VERSION
calculate_mean_abs_smd <- function(W_matrix, X) {
  smd_values <- apply(W_matrix, MARGIN = 2, FUN = function(w) {
    abs(smd(X, w))  # Get absolute SMD values for each covariate
  })
  # smd_values is now a matrix: rows = covariates, columns = replications
  # Take the mean across all covariates and replications
  return(mean(smd_values, na.rm = TRUE))
}

mean_smd_simple <- calculate_mean_abs_smd(W_simple, X)
mean_smd_rer <- calculate_mean_abs_smd(W_rer, X)
mean_smd_ridgerer <- calculate_mean_abs_smd(W_ridgerer, X)
mean_smd_pcarer <- calculate_mean_abs_smd(W_pcarer, X)
mean_smd_pwdrer <- calculate_mean_abs_smd(W_pwdrer, X)
mean_smd_hcm <- calculate_mean_abs_smd(W_hcm, X)
mean_smd_hcm_rer <- calculate_mean_abs_smd(W_hcm_rer, X)

# Create results table with requested columns only
results_table <- data.frame(
  Methods = c("Simple Randomization", "ReR", "RidgeReR", "PCA-ReR", "PWD-ReR", "HCM", "HCM-ReR"),
  mean_bias = c(bias_simple, bias_rer, bias_ridgerer, bias_pcarer, bias_pwdrer, bias_hcm, bias_hcm_rer),
  rmse = c(rmse_simple, rmse_rer, rmse_ridgerer, rmse_pcarer, rmse_pwdrer, rmse_hcm, rmse_hcm_rer),
  mean_smd = c(mean_smd_simple, mean_smd_rer, mean_smd_ridgerer, mean_smd_pcarer, mean_smd_pwdrer, mean_smd_hcm, mean_smd_hcm_rer),
  mean_mahalanobis = c(
    mean(maha_dist_simple, na.rm = TRUE),
    mean(maha_dist_rer, na.rm = TRUE),
    mean(maha_dist_ridgerer, na.rm = TRUE),
    mean(maha_dist_pcarer, na.rm = TRUE),
    mean(maha_dist_pwdrer, na.rm = TRUE),
    mean(maha_dist_hcm, na.rm = TRUE),
    mean(maha_dist_hcm_rer, na.rm = TRUE)
  ),
  mean_iterations = c(
    1, mean(ii_vec_rer), mean(ii_vec_ridgerer), mean(ii_vec_pcarer), 
    mean(ii_vec_pwdrer), 1, mean(ii_vec_hcm_rer)
  )
)

# Print results
cat("\n=== SELECTED METHODS RESULTS ===\n")
print(round(results_table, 4))

# ============================================================================
# VISUALIZATIONS (SELECTED METHODS ONLY)
# ============================================================================

# 1. Treatment Effect Distribution Comparison
tau_df <- data.frame(
  value = c(tau_simple, tau_rer, tau_ridgerer, tau_pcarer, tau_pwdrer, tau_hcm, tau_hcm_rer),
  Methods = factor(rep(c('Simple Randomization', 'ReR', 'RidgeReR', 'PCA-ReR', 'PWD-ReR', 'HCM', 'HCM-ReR'), 
                       rep(reps, 7)),
                   levels = c('Simple Randomization', 'ReR', 'RidgeReR', 'PCA-ReR', 'PWD-ReR', 'HCM', 'HCM-ReR'))
)

tau_density_plot <- ggplot(tau_df, aes(x = value, color = Methods, fill = Methods)) +
  geom_density(alpha = 0.3) +
  geom_vline(aes(xintercept = tau_true), linetype = "dashed", color = "red") +
  xlab(TeX('$\\hat{\\tau}$')) +
  ylab('Density') +
  theme_bw() +
  theme(legend.position = "bottom") +
  annotate("text", x = tau_true, y = Inf, label = "True Effect", 
           vjust = 1.2, color = "red", size = 6)
tau_density_plot <- ggplot(tau_df, aes(x = value, color = Methods, fill = Methods)) +
  geom_density(alpha = 0.3) +
  geom_vline(aes(xintercept = tau_true), linetype = "dashed", color = "red") +
  xlab(TeX('$\\hat{\\tau}$')) +                  
  ylab('Density') +                               
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1, "cm")          
  ) +
  annotate("text", x = tau_true, y = Inf, label = "True Effect", 
           vjust = 1.2, color = "red", size = 6)

tau_density_plot
print(tau_density_plot)

# 2. Performance Comparison
performance_melted <- reshape2::melt(results_table, 
                                     id.vars = "Methods",
                                     variable.name = "Metric",
                                     value.name = "Value")

performance_plot <- ggplot(performance_melted, aes(x = Methods, y = Value, fill = Methods)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Performance Comparison Across Selected Methods",
       x = "Methods", y = "Value")

print(performance_plot)

# 3. Bias vs Balance Trade-off
tradeoff_data <- data.frame(
  Methods = results_table$Methods,
  Bias = abs(results_table$mean_bias),
  Balance = results_table$mean_mahalanobis
)

tradeoff_plot <- ggplot(tradeoff_data, aes(x = Balance, y = Bias, color = Methods)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Methods), vjust = -0.5, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "Bias vs Balance Trade-off",
       subtitle = "Lower left is optimal (low bias, low imbalance)",
       x = "Mean Mahalanobis Distance (Balance)",
       y = "Absolute Mean Bias",
       color = "Methods") +
  theme(legend.position = "bottom")

print(tradeoff_plot)

# ============================================================================
# SAVE RESULTS
# ============================================================================

# Save results table
write.csv(results_table, 
          file = file.path(save_dir_folder, 'selected_methods_results.csv'), 
          row.names = FALSE)
# Save plots
ggsave(file.path(save_dir_folder, 'treatment_effect_distributions.pdf'),
       plot = tau_density_plot, width = 12, height = 8)
ggsave(file.path(save_dir_folder, 'performance_comparison.pdf'),
       plot = performance_plot, width = 14, height = 10)
ggsave(file.path(save_dir_folder, 'bias_balance_tradeoff.pdf'),
       plot = tradeoff_plot, width = 10, height = 8)

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")
cat("Best Balance (lowest Mahalanobis distance):", 
    results_table$Methods[which.min(results_table$mean_mahalanobis)], "\n")
cat("Best Accuracy (lowest RMSE):", 
    results_table$Methods[which.min(results_table$rmse)], "\n")
cat("Most Efficient (fewest iterations):", 
    results_table$Methods[which.min(results_table$mean_iterations)], "\n")
cat("Lowest Bias (closest to 0):", 
    results_table$Methods[which.min(abs(results_table$mean_bias))], "\n")
cat("Best SMD (lowest mean SMD):", 
    results_table$Methods[which.min(results_table$mean_smd)], "\n")

# Ranking by overall performance
results_table$Overall_Score <- 
  (1 - abs(results_table$mean_bias)/max(abs(results_table$mean_bias))) * 0.25 +
  (1 - results_table$rmse/max(results_table$rmse)) * 0.25 +
  (1 - results_table$mean_smd/max(results_table$mean_smd)) * 0.25 +
  (1 - results_table$mean_mahalanobis/max(results_table$mean_mahalanobis)) * 0.25

ranking <- results_table[order(-results_table$Overall_Score), 
                         c("Methods", "Overall_Score")]
cat("\n=== OVERALL RANKING ===\n")
print(ranking)

# HCM-ReR performance highlights
hcm_rer_row <- results_table[results_table$Methods == "HCM-ReR", ]
cat("\n=== HCM-ReR PERFORMANCE HIGHLIGHTS ===\n")
cat("HCM-ReR Rank:", which(ranking$Methods == "HCM-ReR"), "out of", nrow(ranking), "\n")
cat("Bias vs Simple:", round((hcm_rer_row$mean_bias - results_table$mean_bias[1]), 4), "\n")
cat("RMSE vs Simple:", round((hcm_rer_row$rmse - results_table$rmse[1]), 4), "\n")
cat("Balance vs Simple:", round((hcm_rer_row$mean_mahalanobis - results_table$mean_mahalanobis[1]), 2), "\n")

cat("\n=== STUDY COMPLETED SUCCESSFULLY ===\n")
cat("Results saved to:", save_dir_folder, "\n")

# ============================================================================
# COVARIATE BALANCE HEATMAP
# ============================================================================

# Calculate relative variance reduction for each method
r_mdiff_rer <- (apply(cov_diff_simple, 2, var) - apply(cov_diff_rer, 2, var)) / 
  (apply(cov_diff_simple, 2, var) + 1e-30)
r_mdiff_ridgerer <- (apply(cov_diff_simple, 2, var) - apply(cov_diff_ridgerer, 2, var)) / 
  (apply(cov_diff_simple, 2, var) + 1e-30)
r_mdiff_pcarer <- (apply(cov_diff_simple, 2, var) - apply(cov_diff_pcarer, 2, var)) / 
  (apply(cov_diff_simple, 2, var) + 1e-30)
r_mdiff_pwdrer <- (apply(cov_diff_simple, 2, var) - apply(cov_diff_pwdrer, 2, var)) / 
  (apply(cov_diff_simple, 2, var) + 1e-30)
r_mdiff_hcm <- (apply(cov_diff_simple, 2, var) - apply(cov_diff_hcm, 2, var)) / 
  (apply(cov_diff_simple, 2, var) + 1e-30)
r_mdiff_hcm_rer <- (apply(cov_diff_simple, 2, var) - apply(cov_diff_hcm_rer, 2, var)) / 
  (apply(cov_diff_simple, 2, var) + 1e-30)

# Combine methods in desired order
r_mdiff_all <- cbind(r_mdiff_rer, r_mdiff_ridgerer, r_mdiff_pcarer, r_mdiff_pwdrer, r_mdiff_hcm, r_mdiff_hcm_rer)
colnames(r_mdiff_all) <- c("ReR", "Ridge-ReR", "PCA-ReR", "PWD-ReR", "HCM", "HCM-ReR")

# Reorder columns so ReR is first and HCM-ReR is last in the display
r_mdiff_all_reordered <- r_mdiff_all[, rev(colnames(r_mdiff_all))]

# Plot first 25 covariates for visibility
n_covs_to_show <- min(25, nrow(r_mdiff_all_reordered))
r_mdiff_plot <- r_mdiff_all_reordered[1:n_covs_to_show, ]

# Create proper mathematical notation for x-axis labels
x_labels_parse <- paste0("X[", 1:n_covs_to_show, "]")

plt_mdiff_heatmap <- levelplot(r_mdiff_plot,
                               xlab = 'Covariates',
                               ylab = 'Methods',
                               main = 'Covariate Balance Improvement vs Simple Randomization',
                               aspect = 0.4,
                               at = seq(0, 1, length.out = 21),
                               scales = list(
                                 x = list(
                                   at = 1:n_covs_to_show,
                                   labels = parse(text = x_labels_parse),
                                   cex = 0.8
                                 ),
                                 y = list(cex = 0.9)
                               ))
plt_mdiff_heatmap <- levelplot(r_mdiff_plot,
                               xlab = 'Covariates',
                               ylab = 'Methods',
                               aspect = 0.4,
                               at = seq(0, 1, length.out = 21),
                               scales = list(
                                 x = list(
                                   at = 1:n_covs_to_show,
                                   labels = parse(text = x_labels_parse),
                                   cex = 1.3
                                 ),
                                 y = list(
                                   cex = 1.3
                                 )
                               ),
                               colorkey = list(
                                 labels = list(cex = 1.4)
                               ),
                               par.settings = list(
                                 axis.title = list(cex = 1.6),
                                 par.xlab.text = list(cex = 1.6),   
                                 par.ylab.text = list(cex = 1.6),   
                                 key.text = list(cex = 1.6),
                                 key.title = list(cex = 1.4)
                               ))
print(plt_mdiff_heatmap)
