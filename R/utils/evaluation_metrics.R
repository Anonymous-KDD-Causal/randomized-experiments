#' Evaluation Metrics for Randomization Methods
#' 
#' This file contains functions for evaluating the performance of different
#' randomization methods in terms of balance, efficiency, and treatment effect estimation.

#' @import ggplot2
#' @import MASS

#' Calculate Balance Metrics Suite
#' 
#' Computes a comprehensive set of balance metrics for a given allocation.
#' 
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' 
#' @return List containing various balance metrics
#' @export
calculate_balance_metrics <- function(X, w) {
  n <- length(w)
  n1 <- sum(w)
  n0 <- n - n1
  
  if (n1 == 0 || n0 == 0) {
    return(list(
      mahalanobis_distance = Inf,
      max_smd = Inf,
      mean_smd = Inf,
      median_smd = Inf,
      percent_smd_below_01 = 0,
      percent_smd_below_025 = 0,
      max_abs_mean_diff = Inf,
      frobenius_norm = Inf,
      total_variation_distance = 1
    ))
  }
  
  # Mahalanobis distance
  maha_dist_val <- maha_dist(X, w)
  
  # Standardized mean differences
  smd_values <- smd(X, w)
  max_smd_val <- max(abs(smd_values), na.rm = TRUE)
  mean_smd_val <- mean(abs(smd_values), na.rm = TRUE)
  median_smd_val <- median(abs(smd_values), na.rm = TRUE)
  
  # Proportion of SMDs below thresholds
  percent_smd_below_01 <- mean(abs(smd_values) < 0.1, na.rm = TRUE)
  percent_smd_below_025 <- mean(abs(smd_values) < 0.25, na.rm = TRUE)
  
  # Mean differences
  mean_diff_values <- mean_diff(X, w)
  max_abs_mean_diff <- max(abs(mean_diff_values), na.rm = TRUE)
  
  # Frobenius norm of covariance difference
  X_treated <- X[w == 1, , drop = FALSE]
  X_control <- X[w == 0, , drop = FALSE]
  
  if (nrow(X_treated) > 1 && nrow(X_control) > 1) {
    cov_treated <- cov(X_treated)
    cov_control <- cov(X_control)
    frobenius_norm <- norm(cov_treated - cov_control, type = "F")
  } else {
    frobenius_norm <- Inf
  }
  
  # Total variation distance (simplified version)
  total_variation_distance <- 0.5 * sum(abs(colMeans(X_treated) - colMeans(X_control)))
  
  return(list(
    mahalanobis_distance = maha_dist_val,
    max_smd = max_smd_val,
    mean_smd = mean_smd_val,
    median_smd = median_smd_val,
    percent_smd_below_01 = percent_smd_below_01,
    percent_smd_below_025 = percent_smd_below_025,
    max_abs_mean_diff = max_abs_mean_diff,
    frobenius_norm = frobenius_norm,
    total_variation_distance = total_variation_distance
  ))
}

#' Calculate Efficiency Metrics
#' 
#' Computes efficiency-related metrics for randomization methods.
#' 
#' @param result Method result object
#' @param runtime_seconds Runtime in seconds (optional)
#' 
#' @return List containing efficiency metrics
#' @export
calculate_efficiency_metrics <- function(result, runtime_seconds = NULL) {
  n <- length(result$allocation)
  n1 <- sum(result$allocation)
  n0 <- n - n1
  
  # Sample size efficiency (how close to optimal 1:1 allocation)
  sample_size_efficiency <- (4 * n1 * n0) / (n^2)
  
  # Allocation balance
  allocation_balance <- abs(n1 - n0) / n
  
  # Iteration efficiency
  iterations <- ifelse("iterations" %in% names(result), result$iterations, 1)
  iteration_efficiency <- 1 / iterations
  
  # Budget utilization
  budget_utilization <- if ("indicator" %in% names(result) && "ii" %in% names(result)) {
    # This assumes there was a budget and ii represents iterations used
    result$ii / 1000  # Assuming default budget of 1000
  } else {
    iterations / 1000
  }
  
  metrics <- list(
    sample_size_efficiency = sample_size_efficiency,
    allocation_balance = allocation_balance,
    iteration_efficiency = iteration_efficiency,
    budget_utilization = budget_utilization,
    iterations_used = iterations
  )
  
  if (!is.null(runtime_seconds)) {
    metrics$runtime_seconds = runtime_seconds
    metrics$runtime_per_iteration = runtime_seconds / iterations
  }
  
  return(metrics)
}

#' Comprehensive Method Evaluation
#' 
#' Performs comprehensive evaluation of a randomization method result.
#' 
#' @param X Covariate matrix
#' @param result Method result object
#' @param y Outcome vector (optional)
#' @param tau_true True treatment effect (optional)
#' @param runtime_seconds Runtime in seconds (optional)
#' 
#' @return List containing all evaluation metrics
#' @export
comprehensive_evaluation <- function(X, result, y = NULL, tau_true = NULL, runtime_seconds = NULL) {
  # Balance metrics
  balance_metrics <- calculate_balance_metrics(X, result$allocation)
  
  # Efficiency metrics
  efficiency_metrics <- calculate_efficiency_metrics(result, runtime_seconds)
  
  # Treatment effect metrics (if outcomes provided)
  if (!is.null(y)) {
    tau_hat <- tau_est(y, result$allocation)
    
    te_metrics <- list(tau_hat = tau_hat)
    
    if (!is.null(tau_true)) {
      bias_metrics <- calculate_bias(tau_hat, tau_true)
      te_metrics <- c(te_metrics, bias_metrics)
      
      # MSE and RMSE
      te_metrics$mse <- bias_metrics$bias^2
      te_metrics$rmse <- sqrt(te_metrics$mse)
    }
  } else {
    te_metrics <- list()
  }
  
  # Combine all metrics
  all_metrics <- c(balance_metrics, efficiency_metrics, te_metrics)
  
  # Add method information
  all_metrics$method <- ifelse("method" %in% names(result), result$method, "Unknown")
  all_metrics$budget_exhausted <- ifelse("indicator" %in% names(result), result$indicator, FALSE)
  
  return(all_metrics)
}

#' Compare Multiple Methods
#' 
#' Compares multiple randomization methods across various metrics.
#' 
#' @param X Covariate matrix
#' @param results_list List of method results
#' @param method_names Vector of method names
#' @param y Outcome vector (optional)
#' @param tau_true True treatment effect (optional)
#' @param runtime_list List of runtimes (optional)
#' 
#' @return Data frame with comparison results
#' @export
compare_methods <- function(X, results_list, method_names, y = NULL, tau_true = NULL, runtime_list = NULL) {
  if (length(results_list) != length(method_names)) {
    stop("Length of results_list must match length of method_names")
  }
  
  comparison_results <- data.frame()
  
  for (i in seq_along(results_list)) {
    runtime <- if (!is.null(runtime_list)) runtime_list[[i]] else NULL
    
    evaluation <- comprehensive_evaluation(
      X = X, 
      result = results_list[[i]], 
      y = y, 
      tau_true = tau_true, 
      runtime_seconds = runtime
    )
    
    # Convert to data frame row
    eval_df <- data.frame(
      method = method_names[i],
      mahalanobis_distance = evaluation$mahalanobis_distance,
      max_smd = evaluation$max_smd,
      mean_smd = evaluation$mean_smd,
      balance_score = evaluation$percent_smd_below_01,
      sample_size_efficiency = evaluation$sample_size_efficiency,
      iterations = evaluation$iterations_used,
      budget_exhausted = evaluation$budget_exhausted,
      stringsAsFactors = FALSE
    )
    
    # Add treatment effect metrics if available
    if ("tau_hat" %in% names(evaluation)) {
      eval_df$tau_hat <- evaluation$tau_hat
      if ("bias" %in% names(evaluation)) {
        eval_df$bias <- evaluation$bias
        eval_df$abs_bias <- evaluation$abs_bias
        eval_df$rmse <- evaluation$rmse
      }
    }
    
    # Add runtime if available
    if ("runtime_seconds" %in% names(evaluation)) {
      eval_df$runtime_seconds <- evaluation$runtime_seconds
    }
    
    comparison_results <- rbind(comparison_results, eval_df)
  }
  
  return(comparison_results)
}

#' Create Performance Rankings
#' 
#' Creates performance rankings for different methods across multiple criteria.
#' 
#' @param comparison_df Data frame from compare_methods
#' @param criteria List of criteria with weights (e.g., list(balance = 0.4, efficiency = 0.3, bias = 0.3))
#' 
#' @return Data frame with rankings
#' @export
create_performance_rankings <- function(comparison_df, criteria = list(balance = 0.4, efficiency = 0.3, bias = 0.3)) {
  rankings_df <- comparison_df
  
  # Normalize metrics (lower is better for most)
  if ("mahalanobis_distance" %in% colnames(rankings_df)) {
    rankings_df$balance_score_norm <- 1 - (rankings_df$mahalanobis_distance - min(rankings_df$mahalanobis_distance, na.rm = TRUE)) /
      (max(rankings_df$mahalanobis_distance, na.rm = TRUE) - min(rankings_df$mahalanobis_distance, na.rm = TRUE) + 1e-10)
  }
  
  if ("iterations" %in% colnames(rankings_df)) {
    rankings_df$efficiency_score_norm <- 1 - (rankings_df$iterations - min(rankings_df$iterations, na.rm = TRUE)) /
      (max(rankings_df$iterations, na.rm = TRUE) - min(rankings_df$iterations, na.rm = TRUE) + 1e-10)
  }
  
  if ("abs_bias" %in% colnames(rankings_df)) {
    rankings_df$bias_score_norm <- 1 - (rankings_df$abs_bias - min(rankings_df$abs_bias, na.rm = TRUE)) /
      (max(rankings_df$abs_bias, na.rm = TRUE) - min(rankings_df$abs_bias, na.rm = TRUE) + 1e-10)
  }
  
  # Calculate overall score
  rankings_df$overall_score <- 0
  
  if ("balance_score_norm" %in% colnames(rankings_df) && "balance" %in% names(criteria)) {
    rankings_df$overall_score <- rankings_df$overall_score + criteria$balance * rankings_df$balance_score_norm
  }
  
  if ("efficiency_score_norm" %in% colnames(rankings_df) && "efficiency" %in% names(criteria)) {
    rankings_df$overall_score <- rankings_df$overall_score + criteria$efficiency * rankings_df$efficiency_score_norm
  }
  
  if ("bias_score_norm" %in% colnames(rankings_df) && "bias" %in% names(criteria)) {
    rankings_df$overall_score <- rankings_df$overall_score + criteria$bias * rankings_df$bias_score_norm
  }
  
  # Create rankings
  rankings_df$overall_rank <- rank(-rankings_df$overall_score, ties.method = "min")
  
  if ("mahalanobis_distance" %in% colnames(rankings_df)) {
    rankings_df$balance_rank <- rank(rankings_df$mahalanobis_distance, ties.method = "min")
  }
  
  if ("iterations" %in% colnames(rankings_df)) {
    rankings_df$efficiency_rank <- rank(rankings_df$iterations, ties.method = "min")
  }
  
  if ("abs_bias" %in% colnames(rankings_df)) {
    rankings_df$bias_rank <- rank(rankings_df$abs_bias, ties.method = "min")
  }
  
  # Sort by overall rank
  rankings_df <- rankings_df[order(rankings_df$overall_rank), ]
  
  return(rankings_df)
}

#' Visualize Method Comparison
#' 
#' Creates visualization plots for method comparison.
#' 
#' @param comparison_df Data frame from compare_methods
#' @param save_plots Whether to save plots to files
#' @param plot_dir Directory to save plots (if save_plots = TRUE)
#' 
#' @return List of ggplot objects
#' @export
visualize_method_comparison <- function(comparison_df, save_plots = FALSE, plot_dir = "plots") {
  if (!require(ggplot2)) {
    stop("ggplot2 package required for visualization")
  }
  
  if (save_plots && !dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  plots <- list()
  
  # Balance comparison
  if ("mahalanobis_distance" %in% colnames(comparison_df)) {
    p1 <- ggplot(comparison_df, aes(x = reorder(method, mahalanobis_distance), 
                                    y = mahalanobis_distance)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Balance Comparison (Mahalanobis Distance)",
           subtitle = "Lower is better",
           x = "Method", y = "Mahalanobis Distance") +
      scale_y_log10()
    
    plots$balance <- p1
    if (save_plots) ggsave(file.path(plot_dir, "balance_comparison.png"), p1, width = 10, height = 6)
  }
  
  # Efficiency comparison
  if ("iterations" %in% colnames(comparison_df)) {
    p2 <- ggplot(comparison_df, aes(x = reorder(method, iterations), y = iterations)) +
      geom_col(fill = "coral", alpha = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Efficiency Comparison (Iterations)",
           subtitle = "Lower is better",
           x = "Method", y = "Iterations") +
      scale_y_log10()
    
    plots$efficiency <- p2
    if (save_plots) ggsave(file.path(plot_dir, "efficiency_comparison.png"), p2, width = 10, height = 6)
  }
  
  # Bias comparison (if available)
  if ("abs_bias" %in% colnames(comparison_df)) {
    p3 <- ggplot(comparison_df, aes(x = reorder(method, abs_bias), y = abs_bias)) +
      geom_col(fill = "lightgreen", alpha = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Bias Comparison (Absolute Bias)",
           subtitle = "Lower is better",
           x = "Method", y = "Absolute Bias")
    
    plots$bias <- p3
    if (save_plots) ggsave(file.path(plot_dir, "bias_comparison.png"), p3, width = 10, height = 6)
  }
  
  # Combined scatter plot
  if ("mahalanobis_distance" %in% colnames(comparison_df) && "iterations" %in% colnames(comparison_df)) {
    p4 <- ggplot(comparison_df, aes(x = iterations, y = mahalanobis_distance, 
                                    color = method, size = 2)) +
      geom_point(alpha = 0.8) +
      scale_x_log10() +
      scale_y_log10() +
      theme_minimal() +
      labs(title = "Efficiency vs Balance Trade-off",
           subtitle = "Lower left is optimal (low iterations, low Mahalanobis distance)",
           x = "Iterations (log scale)", y = "Mahalanobis Distance (log scale)",
           color = "Method") +
      guides(size = "none")
    
    plots$tradeoff <- p4
    if (save_plots) ggsave(file.path(plot_dir, "efficiency_balance_tradeoff.png"), p4, width = 12, height = 8)
  }
  
  return(plots)
}

#' Statistical Significance Tests
#' 
#' Performs statistical tests to compare method performance.
#' 
#' @param metrics_list List of metric vectors for each method
#' @param method_names Vector of method names
#' @param test_type Type of test ("wilcox" for Wilcoxon, "t" for t-test)
#' 
#' @return Matrix of p-values for pairwise comparisons
#' @export
statistical_significance_tests <- function(metrics_list, method_names, test_type = "wilcox") {
  n_methods <- length(metrics_list)
  p_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods,
                     dimnames = list(method_names, method_names))
  
  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      if (i != j) {
        if (test_type == "wilcox") {
          test_result <- wilcox.test(metrics_list[[i]], metrics_list[[j]])
        } else if (test_type == "t") {
          test_result <- t.test(metrics_list[[i]], metrics_list[[j]])
        }
        p_matrix[i, j] <- test_result$p.value
      } else {
        p_matrix[i, j] <- 1  # Same method
      }
    }
  }
  
  return(p_matrix)
}