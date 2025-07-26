#' Compare Multiple Randomization Methods
#'
#' Convenience function to run and compare HCM, HCM-ReR, and PWD-ReR on the same dataset.
#'
#' @param X Covariate matrix
#' @param methods Vector of method names to compare (default: c("HCM", "HCM-ReR", "PWD-ReR"))
#' @param pa Acceptance probability for re-randomization methods (default: 0.2)
#' @param max_iterations Maximum iterations for re-randomization methods (default: 1000)
#' @param seed Random seed (default: 2020)
#' @param ... Additional arguments passed to individual methods
#'
#' @return List containing results from all methods and a comparison summary
#'
#' @examples
#' X <- generate_covariates(n = 100, p = 10, rho = 0.3, seed = 123)
#' comparison <- compare_methods(X, seed = 123)
#' print(comparison$summary)
#'
#' @export
compare_methods <- function(X, methods = c("HCM", "HCM-ReR", "PWD-ReR"), 
                            pa = 0.2, max_iterations = 1000, seed = 2020, ...) {
  
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

#' Plot Method Comparison
#'
#' Create visualization comparing the performance of different methods.
#'
#' @param comparison_result Result from compare_methods function
#' @param metrics Vector of metrics to plot (default: c("Balance_Score", "Max_SMD", "Mahalanobis_Distance"))
#'
#' @return ggplot object (if ggplot2 is available)
#'
#' @export
plot_method_comparison <- function(comparison_result, metrics = c("Balance_Score", "Max_SMD", "Mahalanobis_Distance")) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }
  
  library(ggplot2)
  
  summary_data <- comparison_result$summary
  
  # Reshape data for plotting
  plot_data <- data.frame()
  for (metric in metrics) {
    if (metric %in% colnames(summary_data)) {
      temp_data <- data.frame(
        Method = summary_data$Method,
        Metric = metric,
        Value = summary_data[[metric]],
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, temp_data)
    }
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Method, y = Value, fill = Method)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    facet_wrap(~Metric, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Method Comparison",
         x = "Method",
         y = "Value") +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  return(p)
}

#' Calculate Pairwise Distances within Matched Pairs
#'
#' For pair-based methods, calculate statistics about within-pair distances.
#'
#' @param result Result from HCM or HCM-ReR
#' @param X Original covariate matrix
#'
#' @return List with pair distance statistics
#'
#' @export
analyze_pair_distances <- function(result, X) {
  
  if (!"matched_pairs" %in% names(result)) {
    stop("Method does not use matched pairs")
  }
  
  distances <- sapply(result$matched_pairs, function(pair) pair$distance)
  
  # Calculate treatment assignment within pairs
  treatment_patterns <- sapply(result$matched_pairs, function(pair) {
    paste(sort(pair$treatments), collapse = "-")
  })
  
  return(list(
    n_pairs = length(distances),
    mean_distance = mean(distances),
    median_distance = median(distances),
    sd_distance = sd(distances),
    min_distance = min(distances),
    max_distance = max(distances),
    treatment_patterns = table(treatment_patterns),
    distance_summary = summary(distances)
  ))
}