#' Calculate Balance Metrics for Treatment Allocation
#'
#' Computes various balance metrics to evaluate the quality of treatment allocation.
#'
#' @param X Covariate matrix (n x p)
#' @param allocation Treatment allocation vector (0/1)
#'
#' @return A list containing various balance metrics
#'
#' @examples
#' X <- generate_covariates(100, 5, seed = 123)
#' allocation <- sample(c(0,1), 100, replace = TRUE)
#' metrics <- calculate_balance_metrics(X, allocation)
#'
#' @export
#' @importFrom MASS ginv
calculate_balance_metrics <- function(X, allocation) {
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
  
  # Mahalanobis distance
  maha_dist_val <- calculate_mahalanobis_distance(X, allocation)
  
  # Standardized mean differences
  smd_values <- calculate_smd(X, allocation)
  max_smd_val <- max(abs(smd_values), na.rm = TRUE)
  mean_smd_val <- mean(abs(smd_values), na.rm = TRUE)
  
  # Balance score (proportion of SMDs < 0.1)
  balance_score <- mean(abs(smd_values) < 0.1, na.rm = TRUE)
  
  # Allocation balance
  allocation_balance <- abs(n1 - n0) / n
  
  # Sample size efficiency
  sample_size_efficiency <- (4 * n1 * n0) / (n^2)
  
  return(list(
    mahalanobis_distance = maha_dist_val,
    max_smd = max_smd_val,
    mean_smd = mean_smd_val,
    balance_score = balance_score,
    allocation_balance = allocation_balance,
    sample_size_efficiency = sample_size_efficiency,
    smd_values = smd_values
  ))
}

#' Calculate Mahalanobis Distance
#'
#' @param X Covariate matrix
#' @param allocation Treatment allocation vector
#' @param lambda Regularization parameter (default: 0)
#'
#' @return Mahalanobis distance value
#'
#' @export
#' @importFrom MASS ginv
calculate_mahalanobis_distance <- function(X, allocation, lambda = 0) {
  n1 <- sum(allocation)
  n0 <- sum(1 - allocation)
  n <- n0 + n1
  
  if (n1 == 0 || n0 == 0) return(Inf)
  
  # Group means
  xbar1 <- t(X) %*% allocation / n1
  xbar0 <- t(X) %*% (1 - allocation) / n0
  delta <- xbar1 - xbar0
  
  # Covariance matrix
  cov_mat <- (n / (n1 * n0 * (n - 1))) * t(X) %*% X
  
  # Regularized inverse
  inv_cov <- MASS::ginv(cov_mat + lambda * diag(nrow(cov_mat)))
  
  # Mahalanobis distance
  mdist <- as.numeric(t(delta) %*% inv_cov %*% delta)
  
  return(mdist)
}

#' Calculate Standardized Mean Differences
#'
#' @param X Covariate matrix
#' @param allocation Treatment allocation vector
#'
#' @return Vector of standardized mean differences
#'
#' @export
calculate_smd <- function(X, allocation) {
  treated <- X[allocation == 1, , drop = FALSE]
  control <- X[allocation == 0, , drop = FALSE]
  
  if (nrow(treated) == 0 || nrow(control) == 0) {
    return(rep(Inf, ncol(X)))
  }
  
  mean_t <- colMeans(treated)
  mean_c <- colMeans(control)
  
  var_t <- apply(treated, 2, var)
  var_c <- apply(control, 2, var)
  pooled_sd <- sqrt((var_t + var_c) / 2)
  
  # Avoid division by zero
  pooled_sd[pooled_sd == 0] <- 1
  
  smd_values <- (mean_t - mean_c) / pooled_sd
  
  return(smd_values)
}

#' Evaluate Treatment Allocation
#'
#' Comprehensive evaluation of a treatment allocation including balance metrics
#' and treatment effect estimation.
#'
#' @param X Covariate matrix
#' @param allocation Treatment allocation vector
#' @param true_tau True treatment effect (for bias calculation)
#' @param beta Coefficient vector for response generation
#' @param sigma Error standard deviation
#' @param response_type Type of response function
#' @param seed Random seed for response generation
#'
#' @return List with evaluation results
#'
#' @export
evaluate_allocation <- function(X, allocation, true_tau = 1, beta = NULL,
                                sigma = 1, response_type = "linear", seed = NULL) {
  
  # Calculate balance metrics
  balance_metrics <- calculate_balance_metrics(X, allocation)
  
  # Generate response and estimate treatment effect
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(beta)) beta <- rep(0.5, ncol(X))
  
  y <- generate_response(X, allocation, beta, true_tau, sigma, response_type, seed)
  
  # Estimate treatment effect
  tau_hat <- mean(y[allocation == 1]) - mean(y[allocation == 0])
  
  # Calculate standard error
  var_treated <- var(y[allocation == 1])
  var_control <- var(y[allocation == 0])
  n_treated <- sum(allocation)
  n_control <- sum(1 - allocation)
  
  if (n_treated > 1 && n_control > 1) {
    se_tau <- sqrt(var_treated/n_treated + var_control/n_control)
  } else {
    se_tau <- Inf
  }
  
  # Combine results
  result <- c(balance_metrics, list(
    tau_hat = tau_hat,
    se_tau = se_tau,
    bias = tau_hat - true_tau,
    relative_bias = (tau_hat - true_tau) / true_tau,
    abs_bias = abs(tau_hat - true_tau),
    rmse = sqrt((tau_hat - true_tau)^2)
  ))
  
  return(result)
}

#' Evaluate Treatment Allocation
#'
#' Comprehensive evaluation of a treatment allocation including balance metrics
#' and treatment effect estimation.
#'
#' @param X Covariate matrix
#' @param allocation Treatment allocation vector
#' @param true_tau True treatment effect (for bias calculation)
#' @param beta Coefficient vector for response generation
#' @param sigma Error standard deviation
#' @param response_type Type of response function
#' @param seed Random seed for response generation
#'
#' @return List with evaluation results
#'
#' @export
compare_methods <- function(X, methods = c("HCM", "HCM-ReR", "PWD-ReR"), 
                            pa = 0.2, max_iterations = 1000, seed = 2020, 
                            true_tau = 1, beta = NULL, sigma = 1, 
                            response_type = "linear", response_seed = NULL, ...) {
  
  set.seed(seed)
  results <- list()
  
  # Run each method
  for (method in methods) {
    cat("Running", method, "...\n")
    
    start_time <- Sys.time()
    
    if (method == "HCM") {
      if (nrow(X) %% 2 == 0) {
        results[[method]] <- hcm(X, seed = seed, ...)
      } else {
        results[[method]] <- list(error = "HCM requires even sample size")
      }
    } else if (method == "HCM-ReR") {
      if (nrow(X) %% 2 == 0) {
        results[[method]] <- hcm_rer(X, pa = pa, max_iterations = max_iterations, seed = seed, ...)
      } else {
        results[[method]] <- list(error = "HCM-ReR requires even sample size")
      }
    } else if (method == "PWD-ReR") {
      results[[method]] <- pwd_rer(X, pa = pa, max_iterations = max_iterations, seed = seed, ...)
    } else {
      warning("Unknown method: ", method)
      next
    }
    
    results[[method]]$runtime <- as.numeric(Sys.time() - start_time)
  }
  
  # Create comparison summary
  summary_data <- data.frame()
  
  for (method in names(results)) {
    result <- results[[method]]
    
    if ("error" %in% names(result)) {
      next
    }
    
    # Calculate balance metrics
    balance_metrics <- calculate_balance_metrics(X, result$allocation)
    
    # Evaluate allocation
    evaluation <- evaluate_allocation(
      X = X, 
      allocation = result$allocation, 
      true_tau = true_tau, 
      beta = beta, 
      sigma = sigma, 
      response_type = response_type, 
      seed = response_seed
    )
    
    # Add evaluation results to summary
    summary_row <- data.frame(
      Method = method,
      Accepted = ifelse("accepted" %in% names(result), result$accepted, TRUE),
      Iterations = ifelse("iterations" %in% names(result), result$iterations, 1),
      Runtime_sec = result$runtime,
      Balance_Score = balance_metrics$balance_score,
      Max_SMD = balance_metrics$max_smd,
      Mean_SMD = balance_metrics$mean_smd,
      Mahalanobis_Distance = balance_metrics$mahalanobis_distance,
      Sample_Size_Efficiency = balance_metrics$sample_size_efficiency,
      Tau_Hat = evaluation$tau_hat,
      Bias = evaluation$bias,
      Relative_Bias = evaluation$relative_bias,
      RMSE = evaluation$rmse,
      stringsAsFactors = FALSE
    )
    
    summary_data <- rbind(summary_data, summary_row)
  }
  
  # Order by balance score (descending)
  summary_data <- summary_data[order(summary_data$Balance_Score, decreasing = TRUE), ]
  
  return(list(
    results = results,
    summary = summary_data,
    best_method = summary_data$Method[1]
  ))
}