#' Comprehensive Simulation Study for Randomization Methods
#' 
#' This script conducts a comprehensive simulation study comparing
#' HCM, HCM-ReR, and PWD-ReR methods across multiple scenarios.

# Load required libraries and functions
source("R/methods/hcm.R")
source("R/methods/hcm_rer.R")
source("R/methods/pwd_rer.R")
source("R/utils/helper_functions.R")
source("R/utils/evaluation_metrics.R")
source("R/utils/data_generation.R")

# Install and load required packages
required_packages <- c("MASS", "ggplot2", "reshape2", "dplyr", "gridExtra")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#' Run Single Simulation Replication
#' 
#' Runs all methods on a single simulated dataset
#' 
#' @param X Covariate matrix
#' @param outcome_function Function to generate outcomes given allocation
#' @param true_ate True average treatment effect
#' @param simulation_params List of simulation parameters
#' @param methods_to_run Vector of method names to run
#' 
#' @return Data frame with results for this replication
run_single_replication <- function(X, outcome_function, true_ate, 
                                   simulation_params, methods_to_run = "all") {
  
  n <- nrow(X)
  
  # Define which methods to run
  if (methods_to_run == "all") {
    methods_to_run <- c("Simple", "HCM", "HCM-ReR", "PWD-ReR")
  }
  
  results_list <- list()
  runtimes <- list()
  
  # Simple Randomization (baseline)
  if ("Simple" %in% methods_to_run) {
    start_time <- Sys.time()
    simple_result <- simple_randomization(n, seed = simulation_params$rep_seed)
    end_time <- Sys.time()
    
    results_list$Simple <- simple_result
    runtimes$Simple <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  # Hierarchical Clustering Matching (HCM)
  if ("HCM" %in% methods_to_run) {
    start_time <- Sys.time()
    hcm_result <- HCM(X, seed = simulation_params$rep_seed)
    end_time <- Sys.time()
    
    results_list$HCM <- hcm_result
    runtimes$HCM <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  # HCM with Rerandomization (HCM-ReR)
  if ("HCM-ReR" %in% methods_to_run) {
    start_time <- Sys.time()
    hcm_rer_result <- HCM_ReR(
      pa = simulation_params$pa, 
      X = X, 
      var_explained = simulation_params$var_explained,
      n_budget = simulation_params$n_budget, 
      seed = simulation_params$rep_seed
    )
    end_time <- Sys.time()
    
    results_list$`HCM-ReR` <- hcm_rer_result
    runtimes$`HCM-ReR` <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  # PCA-Weighted Distance Rerandomization (PWD-ReR)
  if ("PWD-ReR" %in% methods_to_run) {
    start_time <- Sys.time()
    pwd_rer_result <- PWD_ReR(
      pa = simulation_params$pa, 
      X = X, 
      var_explained = simulation_params$var_explained,
      n_budget = simulation_params$n_budget, 
      seed = simulation_params$rep_seed
    )
    end_time <- Sys.time()
    
    results_list$`PWD-ReR` <- pwd_rer_result
    runtimes$`PWD-ReR` <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  # Generate outcomes for each method and evaluate
  method_names <- names(results_list)
  evaluation_results <- data.frame()
  
  for (method_name in method_names) {
    result <- results_list[[method_name]]
    runtime <- runtimes[[method_name]]
    
    # Generate outcomes for this allocation
    y <- outcome_function(result$allocation)
    
    # Comprehensive evaluation
    evaluation <- comprehensive_evaluation(
      X = X, 
      result = result, 
      y = y, 
      tau_true = true_ate, 
      runtime_seconds = runtime
    )
    
    # Convert to data frame row
    eval_row <- data.frame(
      replication = simulation_params$rep_id,
      method = method_name,
      sample_size = n,
      n_covariates = ncol(X),
      correlation_structure = simulation_params$correlation_structure,
      rho = simulation_params$rho,
      pa = simulation_params$pa,
      true_ate = true_ate,
      
      # Balance metrics
      mahalanobis_distance = evaluation$mahalanobis_distance,
      max_smd = evaluation$max_smd,
      mean_smd = evaluation$mean_smd,
      balance_score = evaluation$percent_smd_below_01,
      
      # Efficiency metrics
      iterations = evaluation$iterations_used,
      runtime_seconds = evaluation$runtime_seconds,
      budget_exhausted = evaluation$budget_exhausted,
      sample_size_efficiency = evaluation$sample_size_efficiency,
      
      # Treatment effect metrics
      tau_hat = evaluation$tau_hat,
      bias = evaluation$bias,
      abs_bias = evaluation$abs_bias,
      relative_bias = evaluation$relative_bias,
      mse = evaluation$mse,
      rmse = evaluation$rmse,
      
      stringsAsFactors = FALSE
    )
    
    evaluation_results <- rbind(evaluation_results, eval_row)
  }
  
  return(evaluation_results)
}

#' Run Simulation Study Configuration
#' 
#' Runs a complete simulation study for a given configuration
#' 
#' @param config Simulation configuration list
#' @param n_replications Number of simulation replications
#' @param save_individual_results Whether to save individual replication results
#' @param results_dir Directory to save results
#' 
#' @return Aggregated simulation results
run_simulation_configuration <- function(config, n_replications = 500, 
                                         save_individual_results = FALSE,
                                         results_dir = "results") {
  
  cat("Running simulation configuration:\n")
  cat("  Sample size:", config$n, "\n")
  cat("  Covariates:", config$p, "\n") 
  cat("  Correlation:", config$correlation_structure, "(rho =", config$rho, ")\n")
  cat("  Effect type:", config$effect_type, "\n")
  cat("  Replications:", n_replications, "\n\n")
  
  # Create results directory if needed
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  all_results <- data.frame()
  
  # Progress tracking
  start_time <- Sys.time()
  
  for (rep in 1:n_replications) {
    # Show progress
    if (rep %% 50 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat("Completed", rep, "/", n_replications, "replications in", 
          round(elapsed, 2), "minutes\n")
    }
    
    # Generate data for this replication
    sim_data <- generate_simulation_data(
      n = config$n,
      p = config$p,
      correlation_structure = config$correlation_structure,
      rho = config$rho,
      effect_type = config$effect_type,
      base_effect = config$base_effect,
      noise_sd = config$noise_sd,
      seed = config$base_seed + rep
    )
    
    # Set up parameters for this replication
    rep_params <- list(
      rep_id = rep,
      rep_seed = config$base_seed + rep,
      pa = config$pa,
      var_explained = config$var_explained,
      n_budget = config$n_budget,
      correlation_structure = config$correlation_structure,
      rho = config$rho
    )
    
    # Run single replication
    tryCatch({
      rep_results <- run_single_replication(
        X = sim_data$X,
        outcome_function = sim_data$generate_outcomes_function,
        true_ate = sim_data$true_ate,
        simulation_params = rep_params
      )
      
      all_results <- rbind(all_results, rep_results)
      
    }, error = function(e) {
      cat("Error in replication", rep, ":", e$message, "\n")
    })
  }
  
  # Save individual results if requested
  if (save_individual_results) {
    config_name <- paste0("sim_n", config$n, "_p", config$p, "_", 
                          config$correlation_structure, "_rho", config$rho)
    write.csv(all_results, 
              file.path(results_dir, paste0(config_name, "_individual_results.csv")),
              row.names = FALSE)
  }
  
  # Calculate summary statistics
  summary_results <- all_results %>%
    group_by(method) %>%
    summarise(
      n_replications = n(),
      
      # Balance metrics
      mean_mahalanobis = mean(mahalanobis_distance, na.rm = TRUE),
      median_mahalanobis = median(mahalanobis_distance, na.rm = TRUE),
      mean_max_smd = mean(max_smd, na.rm = TRUE),
      mean_balance_score = mean(balance_score, na.rm = TRUE),
      
      # Efficiency metrics  
      mean_iterations = mean(iterations, na.rm = TRUE),
      median_iterations = median(iterations, na.rm = TRUE),
      mean_runtime = mean(runtime_seconds, na.rm = TRUE),
      budget_exhausted_rate = mean(budget_exhausted, na.rm = TRUE),
      
      # Treatment effect metrics
      mean_tau_hat = mean(tau_hat, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_abs_bias = mean(abs_bias, na.rm = TRUE),
      median_abs_bias = median(abs_bias, na.rm = TRUE),
      rmse = sqrt(mean(mse, na.rm = TRUE)),
      bias_sd = sd(bias, na.rm = TRUE),
      
      # Configuration info
      sample_size = first(sample_size),
      n_covariates = first(n_covariates),
      correlation_structure = first(correlation_structure),
      rho = first(rho),
      true_ate = first(true_ate),
      
      .groups = 'drop'
    )
  
  # Add ranking information
  summary_results$balance_rank <- rank(summary_results$mean_mahalanobis)
  summary_results$efficiency_rank <- rank(summary_results$mean_iterations)  
  summary_results$bias_rank <- rank(summary_results$mean_abs_bias)
  summary_results$overall_rank <- rank(summary_results$mean_abs_bias + 
                                         0.5 * rank(summary_results$mean_mahalanobis) / nrow(summary_results) +
                                         0.3 * rank(summary_results$mean_iterations) / nrow(summary_results))
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat("Configuration completed in", round(total_time, 2), "minutes\n\n")
  
  return(list(
    summary = summary_results,
    individual_results = all_results,
    configuration = config,
    total_runtime_minutes = total_time
  ))
}

#' Run Complete Simulation Study
#' 
#' Runs the complete simulation study across multiple configurations
#' 
#' @param n_replications Number of replications per configuration
#' @param save_results Whether to save results to files
#' @param results_dir Directory for saving results
#' 
#' @return List containing all simulation results
run_complete_simulation_study <- function(n_replications = 500, 
                                          save_results = TRUE,
                                          results_dir = "simulation_results") {
  
  cat("=== COMPREHENSIVE SIMULATION STUDY ===\n")
  cat("Comparing HCM, HCM-ReR, and PWD-ReR methods\n")
  cat("Replications per configuration:", n_replications, "\n\n")
  
  # Define simulation configurations
  configurations <- list(
    # Small sample, low correlation
    list(
      name = "small_low_corr",
      n = 100, p = 10, correlation_structure = "exchangeable", rho = 0.3,
      effect_type = "linear", base_effect = 1.0, noise_sd = 1.0,
      pa = 0.05, var_explained = 0.7, n_budget = 1000, base_seed = 2020
    ),
    
    # Medium sample, medium correlation  
    list(
      name = "medium_med_corr",
      n = 200, p = 20, correlation_structure = "ar1", rho = 0.5,
      effect_type = "linear", base_effect = 1.0, noise_sd = 1.0,
      pa = 0.05, var_explained = 0.7, n_budget = 1000, base_seed = 3020
    ),
    
    # Large sample, high correlation
    list(
      name = "large_high_corr",
      n = 400, p = 30, correlation_structure = "ar1", rho = 0.8,
      effect_type = "linear", base_effect = 1.0, noise_sd = 1.0,
      pa = 0.05, var_explained = 0.7, n_budget = 1000, base_seed = 4020
    ),
    
    # High dimensional
    list(
      name = "high_dimensional",
      n = 200, p = 50, correlation_structure = "block", rho = c(0.7, 0.5, 0.3),
      effect_type = "linear", base_effect = 1.0, noise_sd = 1.0,
      pa = 0.05, var_explained = 0.7, n_budget = 1000, base_seed = 5020
    ),
    
    # Nonlinear effects
    list(
      name = "nonlinear_effects",
      n = 200, p = 20, correlation_structure = "ar1", rho = 0.5,
      effect_type = "nonlinear", base_effect = 1.0, noise_sd = 1.0,
      pa = 0.05, var_explained = 0.7, n_budget = 1000, base_seed = 6020
    )
  )
  
  # Create results directory
  if (save_results && !dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Run each configuration
  all_configurations_results <- list()
  study_start_time <- Sys.time()
  
  for (i in seq_along(configurations)) {
    config <- configurations[[i]]
    cat("=== Configuration", i, "of", length(configurations), ":", config$name, "===\n")
    
    config_results <- run_simulation_configuration(
      config = config,
      n_replications = n_replications,
      save_individual_results = save_results,
      results_dir = results_dir
    )
    
    all_configurations_results[[config$name]] <- config_results
    
    # Save configuration summary
    if (save_results) {
      write.csv(config_results$summary,
                file.path(results_dir, paste0(config$name, "_summary.csv")),
                row.names = FALSE)
    }
  }
  
  # Create overall summary across configurations
  overall_summary <- data.frame()
  
  for (config_name in names(all_configurations_results)) {
    config_summary <- all_configurations_results[[config_name]]$summary
    config_summary$configuration <- config_name
    overall_summary <- rbind(overall_summary, config_summary)
  }
  
  # Calculate cross-configuration rankings
  method_rankings <- overall_summary %>%
    group_by(method) %>%
    summarise(
      avg_balance_rank = mean(balance_rank),
      avg_efficiency_rank = mean(efficiency_rank),
      avg_bias_rank = mean(bias_rank),
      avg_overall_rank = mean(overall_rank),
      
      avg_mahalanobis = mean(mean_mahalanobis),
      avg_iterations = mean(mean_iterations),
      avg_abs_bias = mean(mean_abs_bias),
      avg_rmse = mean(rmse),
      
      n_configurations = n(),
      .groups = 'drop'
    ) %>%
    arrange(avg_overall_rank)
  
  total_study_time <- as.numeric(difftime(Sys.time(), study_start_time, units = "hours"))
  
  # Save overall results
  if (save_results) {
    write.csv(overall_summary, 
              file.path(results_dir, "overall_summary.csv"), 
              row.names = FALSE)
    write.csv(method_rankings,
              file.path(results_dir, "method_rankings.csv"),
              row.names = FALSE)
  }
  
  # Print final summary
  cat("\n=== SIMULATION STUDY COMPLETED ===\n")
  cat("Total study time:", round(total_study_time, 2), "hours\n")
  cat("Configurations tested:", length(configurations), "\n")
  cat("Total replications:", n_replications * length(configurations), "\n\n")
  
  cat("=== METHOD RANKINGS (Lower is better) ===\n")
  print(method_rankings[, c("method", "avg_overall_rank", "avg_balance_rank", 
                            "avg_efficiency_rank", "avg_bias_rank")])
  
  return(list(
    configuration_results = all_configurations_results,
    overall_summary = overall_summary,
    method_rankings = method_rankings,
    total_runtime_hours = total_study_time
  ))
}

#' Generate Simulation Report
#' 
#' Generates a comprehensive report of simulation results
#' 
#' @param simulation_results Results from run_complete_simulation_study
#' @param save_plots Whether to save plots
#' @param results_dir Directory for saving results
#' 
#' @export
generate_simulation_report <- function(simulation_results, save_plots = TRUE, 
                                       results_dir = "simulation_results") {
  
  if (!require(ggplot2)) stop("ggplot2 required for report generation")
  
  overall_summary <- simulation_results$overall_summary
  method_rankings <- simulation_results$method_rankings
  
  # Create plots directory
  plots_dir <- file.path(results_dir, "plots")
  if (save_plots && !dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Plot 1: Balance comparison across configurations
  p1 <- ggplot(overall_summary, aes(x = configuration, y = mean_mahalanobis, fill = method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Balance Performance Across Configurations",
         subtitle = "Mahalanobis Distance (Lower is Better)",
         x = "Configuration", y = "Mean Mahalanobis Distance (log scale)", fill = "Method")
  
  if (save_plots) ggsave(file.path(plots_dir, "balance_comparison.png"), p1, width = 12, height = 8)
  print(p1)
  
  # Plot 2: Efficiency comparison  
  p2 <- ggplot(overall_summary, aes(x = configuration, y = mean_iterations, fill = method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Efficiency Comparison Across Configurations", 
         subtitle = "Iterations Required (Lower is Better)",
         x = "Configuration", y = "Mean Iterations (log scale)", fill = "Method")
  
  if (save_plots) ggsave(file.path(plots_dir, "efficiency_comparison.png"), p2, width = 12, height = 8)
  print(p2)
  
  # Plot 3: Bias comparison
  p3 <- ggplot(overall_summary, aes(x = configuration, y = mean_abs_bias, fill = method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Treatment Effect Estimation Bias",
         subtitle = "Mean Absolute Bias (Lower is Better)", 
         x = "Configuration", y = "Mean Absolute Bias", fill = "Method")
  
  if (save_plots) ggsave(file.path(plots_dir, "bias_comparison.png"), p3, width = 12, height = 8)
  print(p3)
  
  # Plot 4: Overall performance radar chart (simplified as scatter)
  p4 <- ggplot(method_rankings, aes(x = avg_balance_rank, y = avg_bias_rank, 
                                    color = method, size = 1/avg_overall_rank)) +
    geom_point(alpha = 0.8) +
    geom_text(aes(label = method), vjust = -0.5, hjust = 0.5) +
    theme_minimal() +
    labs(title = "Overall Performance Trade-offs",
         subtitle = "Lower ranks are better (closer to origin is optimal)",
         x = "Average Balance Rank", y = "Average Bias Rank", 
         color = "Method", size = "Overall Performance") +
    guides(size = "none")
  
  if (save_plots) ggsave(file.path(plots_dir, "overall_performance.png"), p4, width = 10, height = 8)
  print(p4)
  
  # Generate text summary
  cat("\n=== SIMULATION STUDY REPORT ===\n\n")
  
  cat("Best Overall Method:", method_rankings$method[1], "\n")
  cat("Best Balance Method:", method_rankings$method[which.min(method_rankings$avg_balance_rank)], "\n")
  cat("Most Efficient Method:", method_rankings$method[which.min(method_rankings$avg_efficiency_rank)], "\n")
  cat("Least Biased Method:", method_rankings$method[which.min(method_rankings$avg_bias_rank)], "\n\n")
  
  cat("=== DETAILED RANKINGS ===\n")
  rankings_print <- method_rankings[, c("method", "avg_overall_rank", "avg_mahalanobis", 
                                        "avg_iterations", "avg_abs_bias", "avg_rmse")]
  print(round(rankings_print, 4))
  
  cat("\n=== RECOMMENDATIONS ===\n")
  if (method_rankings$method[1] == "HCM-ReR") {
    cat("HCM-ReR shows the best overall performance, combining the clustering benefits\n")
    cat("of HCM with the balance optimization of rerandomization.\n")
  } else if (method_rankings$method[1] == "PWD-ReR") {
    cat("PWD-ReR shows the best overall performance, effectively using PCA weights\n") 
    cat("to focus balance optimization on the most important covariate dimensions.\n")
  } else if (method_rankings$method[1] == "HCM") {
    cat("HCM shows the best overall performance through simple but effective\n")
    cat("hierarchical clustering-based matching.\n")
  }
  
  invisible(list(plots = list(p1, p2, p3, p4), rankings = method_rankings))
}

# Main execution function
main_simulation_study <- function() {
  cat("Starting comprehensive simulation study...\n")
  
  # Run the complete study
  results <- run_complete_simulation_study(
    n_replications = 500,  # Adjust based on computational resources
    save_results = TRUE,
    results_dir = "simulation_results"
  )
  
  # Generate report
  generate_simulation_report(
    simulation_results = results,
    save_plots = TRUE,
    results_dir = "simulation_results"
  )
  
  cat("Simulation study completed! Check 'simulation_results' directory for outputs.\n")
  return(results)
}

# Uncomment to run the complete study
# results <- main_simulation_study()