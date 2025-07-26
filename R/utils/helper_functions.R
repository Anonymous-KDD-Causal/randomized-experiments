#' Helper Functions for Randomization Methods
#' 
#' This file contains utility functions used across different randomization methods.

# Load required packages
#' @import MASS
#' @import stats

#' Calculate Mahalanobis Distance
#' 
#' Computes the Mahalanobis distance for treatment allocation balance assessment.
#' 
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' @param lambda Ridge regularization parameter (default: 0)
#' 
#' @return Mahalanobis distance value
#' @export
maha_dist <- function(X, w, lambda = 0) {
  n1 <- sum(w)
  n0 <- sum(1-w)
  n <- n0 + n1
  
  if (n1 == 0 || n0 == 0) {
    return(Inf)
  }
  
  xbar1 <- t(X) %*% w / n1
  xbar0 <- t(X) %*% (1-w) / n0
  delta <- xbar1 - xbar0
  cov_mat <- (n/(n1*n0*(n-1))) * t(X) %*% X
  
  # Use generalized inverse to handle potential singularity
  mdist <- as.numeric(t(delta) %*% MASS::ginv(cov_mat + lambda*diag(nrow(cov_mat))) %*% delta)
  return(mdist)
}

#' Calculate Standardized Mean Differences (SMD)
#' 
#' Computes standardized mean differences for each covariate between treatment groups.
#' 
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' 
#' @return Vector of standardized mean differences
#' @export
smd <- function(X, w) {
  treated <- X[w == 1, , drop = FALSE]
  control <- X[w == 0, , drop = FALSE]
  
  if (nrow(treated) == 0 || nrow(control) == 0) {
    return(rep(Inf, ncol(X)))
  }
  
  mean_t <- colMeans(treated)
  mean_c <- colMeans(control)
  
  var_t <- apply(treated, 2, var)
  var_c <- apply(control, 2, var)
  pooled_sd <- sqrt((var_t + var_c) / 2)
  
  # Handle zero variance case
  pooled_sd[pooled_sd == 0] <- 1
  
  smd_values <- (mean_t - mean_c) / pooled_sd
  return(smd_values)
}

#' Calculate Mean Differences
#' 
#' Computes mean differences for each covariate between treatment groups.
#' 
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' 
#' @return Vector of mean differences
#' @export
mean_diff <- function(X, w) {
  n1 <- sum(w)
  n0 <- sum(1-w)
  
  if (n1 == 0 || n0 == 0) {
    return(rep(Inf, ncol(X)))
  }
  
  xbar1 <- t(X) %*% w / n1
  xbar0 <- t(X) %*% (1-w) / n0
  delta <- xbar1 - xbar0
  return(as.vector(delta))
}

#' Estimate Treatment Effect
#' 
#' Estimates the average treatment effect from outcomes and treatment allocation.
#' 
#' @param y Outcome vector
#' @param w Treatment allocation vector (0/1)
#' 
#' @return Estimated treatment effect
#' @export
tau_est <- function(y, w) {
  if (sum(w) == 0 || sum(1-w) == 0) {
    return(NA)
  }
  mean(y[w==1]) - mean(y[w==0])
}

#' Calculate Treatment Effect Bias
#' 
#' Calculates bias, absolute bias, and relative bias of treatment effect estimates.
#' 
#' @param tau_hat Estimated treatment effect
#' @param tau_true True treatment effect
#' 
#' @return List containing bias, absolute bias, and relative bias
#' @export
calculate_bias <- function(tau_hat, tau_true) {
  bias <- tau_hat - tau_true
  abs_bias <- abs(bias)
  relative_bias <- ifelse(tau_true != 0, bias / tau_true, NA)
  
  return(list(
    bias = bias,
    abs_bias = abs_bias,
    relative_bias = relative_bias
  ))
}

#' Comprehensive Evaluation Metrics
#' 
#' Computes comprehensive evaluation metrics for a given allocation.
#' 
#' @param X Covariate matrix (n x p)
#' @param allocation Treatment allocation vector (0/1)
#' @param tau_true True treatment effect (optional)
#' @param y Outcomes (optional, for treatment effect estimation)
#' 
#' @return List of evaluation metrics
#' @export
evaluate_allocation <- function(X, allocation, tau_true = NULL, y = NULL) {
  n <- length(allocation)
  n1 <- sum(allocation)
  n0 <- n - n1
  
  if (n1 == 0 || n0 == 0) {
    return(list(
      mahalanobis_distance = Inf,
      max_smd = Inf,
      mean_smd = Inf,
      balance_score = 0,
      allocation_balance = abs(n1 - n0) / n,
      sample_size_efficiency = 0
    ))
  }
  
  # Balance metrics
  maha_dist_val <- maha_dist(X, allocation)
  smd_values <- smd(X, allocation)
  max_smd_val <- max(abs(smd_values), na.rm = TRUE)
  mean_smd_val <- mean(abs(smd_values), na.rm = TRUE)
  balance_score <- mean(abs(smd_values) < 0.1, na.rm = TRUE)
  allocation_balance <- abs(n1 - n0) / n
  sample_size_efficiency <- (4 * n1 * n0) / (n^2)
  
  result <- list(
    mahalanobis_distance = maha_dist_val,
    max_smd = max_smd_val,
    mean_smd = mean_smd_val,
    balance_score = balance_score,
    allocation_balance = allocation_balance,
    sample_size_efficiency = sample_size_efficiency
  )
  
  # Treatment effect metrics (if outcomes provided)
  if (!is.null(y)) {
    tau_hat <- tau_est(y, allocation)
    result$tau_hat <- tau_hat
    
    if (!is.null(tau_true)) {
      bias_metrics <- calculate_bias(tau_hat, tau_true)
      result <- c(result, bias_metrics)
    }
  }
  
  return(result)
}

#' Create Comparison Table
#' 
#' Creates a formatted comparison table for multiple allocation methods.
#' 
#' @param results_list List of results from different methods
#' @param method_names Names of the methods
#' 
#' @return Data frame with comparison metrics
#' @export
create_comparison_table <- function(results_list, method_names) {
  if (length(results_list) != length(method_names)) {
    stop("Length of results_list must match length of method_names")
  }
  
  # Extract metrics from each result
  metrics_df <- data.frame(
    Method = method_names,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    
    # Handle different result structures
    if ("allocation" %in% names(result)) {
      # This is a method result with allocation
      metrics_df$Iterations[i] <- ifelse("iterations" %in% names(result), result$iterations, 1)
      metrics_df$Budget_Exhausted[i] <- ifelse("indicator" %in% names(result), result$indicator, FALSE)
      
      if ("matched_pairs" %in% names(result)) {
        metrics_df$Matched_Pairs[i] <- result$n_matched_pairs
        metrics_df$Matching_Rate[i] <- result$n_matched_units / length(result$allocation)
      }
    } else {
      # This is evaluation metrics
      for (metric in names(result)) {
        if (metric %in% colnames(metrics_df)) {
          metrics_df[i, metric] <- result[[metric]]
        } else {
          metrics_df[i, metric] <- result[[metric]]
        }
      }
    }
  }
  
  return(metrics_df)
}

#' Simple Random Allocation
#' 
#' Generates a simple random allocation for comparison purposes.
#' 
#' @param n Sample size
#' @param seed Random seed
#' 
#' @return Treatment allocation vector
#' @export
simple_randomization <- function(n, seed = 2020) {
  set.seed(seed)
  allocation <- sample(c(rep(0, n/2), rep(1, n/2)), n, replace = FALSE)
  
  return(list(
    allocation = allocation,
    method = "Simple",
    iterations = 1
  ))
}

#' Check Input Validity
#' 
#' Validates inputs for randomization methods.
#' 
#' @param X Covariate matrix
#' @param pa Acceptance probability (optional)
#' @param n_budget Budget (optional)
#' 
#' @return TRUE if inputs are valid, stops with error otherwise
#' @export
check_inputs <- function(X, pa = NULL, n_budget = NULL) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame")
  }
  
  if (nrow(X) < 2) {
    stop("X must have at least 2 rows")
  }
  
  if (ncol(X) < 1) {
    stop("X must have at least 1 column")
  }
  
  if (!is.null(pa)) {
    if (pa <= 0 || pa >= 1) {
      stop("pa must be between 0 and 1")
    }
  }
  
  if (!is.null(n_budget)) {
    if (n_budget <= 0 || !is.numeric(n_budget)) {
      stop("n_budget must be a positive number")
    }
  }
  
  # Check for missing values
  if (any(is.na(X))) {
    warning("X contains missing values. Consider imputation or removal.")
  }
  
  return(TRUE)
}

#' Print Method Summary
#' 
#' Prints a formatted summary of method results.
#' 
#' @param result Method result object
#' 
#' @export
print_method_summary <- function(result) {
  cat("=== Method Summary ===\n")
  cat("Method:", ifelse("method" %in% names(result), result$method, "Unknown"), "\n")
  
  if ("allocation" %in% names(result)) {
    cat("Sample size:", length(result$allocation), "\n")
    cat("Treatment group size:", sum(result$allocation), "\n")
    cat("Control group size:", sum(1 - result$allocation), "\n")
  }
  
  if ("iterations" %in% names(result)) {
    cat("Iterations used:", result$iterations, "\n")
  }
  
  if ("indicator" %in% names(result)) {
    cat("Budget exhausted:", result$indicator, "\n")
  }
  
  if ("n_matched_pairs" %in% names(result)) {
    cat("Matched pairs:", result$n_matched_pairs, "\n")
    cat("Matching rate:", round(result$n_matched_units / length(result$allocation), 3), "\n")
  }
  
  if ("k" %in% names(result)) {
    cat("PCA components used:", result$k, "\n")
  }
  
  cat("=====================\n")
}