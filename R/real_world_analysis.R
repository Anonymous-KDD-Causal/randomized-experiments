#' Real-World Data Analysis for Randomization Methods
#' 
#' This script applies HCM, HCM-ReR, and PWD-ReR methods to the IHDP dataset
#' and evaluates their performance on real-world data.

# Load required functions and libraries
source("R/methods/hcm.R")
source("R/methods/hcm_rer.R")
source("R/methods/pwd_rer.R")
source("R/utils/helper_functions.R")
source("R/utils/evaluation_metrics.R")

# Load required packages
required_packages <- c("MASS", "ggplot2", "dplyr", "glmnet", "reshape2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#' Load and Preprocess IHDP Data
#' 
#' Loads the IHDP dataset and performs necessary preprocessing
#' 
#' @param data_path Path to the IHDP data files
#' 
#' @return List containing processed data and metadata
load_and_preprocess_ihdp <- function(data_path = "data/") {
  
  cat("Loading IHDP data...\n")
  
  # Note: This function assumes the IHDP data files are available
  # Users should place the data files in the specified directory
  
  # Load the data files (adjust paths as needed)
  tryCatch({
    load(file.path(data_path, "example.data"))
    load(file.path(data_path, "sim.data"))
  }, error = function(e) {
    stop("IHDP data files not found. Please ensure data files are in the correct directory.\n",
         "Expected files: example.data, sim.data\n",
         "Error: ", e$message)
  })
  
  # Combine data with treatment effect information
  imp1 <- data.frame(iqsb.36 = ihdp$iqsb.36, imp1)
  
  # Define covariate names
  covs.cont.n <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage")
  covs.cat.n <- c("sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                  "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                  "ark", "ein", "har", "mia", "pen", "tex", "was")
  
  p <- length(c(covs.cont.n, covs.cat.n))
  
  # Select relevant variables
  imp1 <- imp1[, c('iqsb.36', 'treat', covs.cont.n, covs.cat.n)]
  
  # Clean categorical variables
  for (cov_name in covs.cat.n) {
    cat(cov_name, "unique values:", unique(imp1[, cov_name]), '\n')
  }
  
  # Fix specific variable encoding
  imp1[, 'first'] <- ifelse(imp1[, 'first'] == 2, 1, 0)
  
  # Remove missing values
  imp1 <- na.omit(imp1)
  ihdp_clean <- imp1
  
  # Rename columns for easier handling
  colnames(ihdp_clean) <- c('y', 'treatment', paste('x', 1:(ncol(ihdp_clean)-2), sep=''))
  
  # Standardize all variables
  ihdp_clean[, 1:27] <- scale(ihdp_clean[, 1:27])
  
  cat("Data preprocessing completed.\n")
  cat("Sample size:", nrow(ihdp_clean), "\n")
  cat("Number of covariates:", ncol(ihdp_clean) - 2, "\n")
  
  return(list(
    data = ihdp_clean,
    continuous_vars = covs.cont.n,
    categorical_vars = covs.cat.n,
    variable_descriptions = list(
      continuous = "Birth weight, head circumference, preterm, birth order, neonatal health, mother age",
      categorical = "Sex, twin, mother married, education levels, smoking, firstborn, alcohol, drugs, work duration, prenatal care, site indicators"
    )
  ))
}

#' Fit Outcome Model for IHDP Data
#' 
#' Fits a regularized linear model to estimate treatment effects
#' 
#' @param ihdp_data Preprocessed IHDP data
#' @param use_interactions Whether to include interaction terms
#' 
#' @return List containing fitted model and estimated treatment effect
fit_ihdp_outcome_model <- function(ihdp_data, use_interactions = TRUE) {
  
  cat("Fitting outcome model...\n")
  
  ihdp <- ihdp_data$data
  n <- nrow(ihdp)
  
  # Create model matrix
  if (use_interactions) {
    # Include all variables and their pairwise interactions
    vars <- paste(paste('x', 1:(ncol(ihdp)-2), sep=''), collapse = '+')
    fml <- as.formula(paste('~ (', vars, ') ^2 - 1'))
  } else {
    # Linear terms only
    vars <- paste(paste('x', 1:(ncol(ihdp)-2), sep=''), collapse = '+')
    fml <- as.formula(paste('~', vars, '- 1'))
  }
  
  design_mat <- data.frame(treatment = ihdp$treatment,
                           model.matrix(fml, data = ihdp))
  design_mat <- scale(design_mat)
  
  # Fit regularized model
  cv_output <- cv.glmnet(x = as.matrix(design_mat), y = ihdp$y, family = 'gaussian')
  best_lam <- cv_output$lambda.min
  lm_fit <- glmnet(x = as.matrix(design_mat), y = ihdp$y, 
                   lambda = best_lam, family = 'gaussian')
  
  # Extract treatment effect coefficient
  tau <- coef(lm_fit)[2]  # Treatment coefficient
  
  cat("Estimated treatment effect (tau):", round(tau, 4), "\n")
  
  # Create outcome generation function
  ihdp_y_gen <- function(X, w) {
    # Combine treatment allocation with covariates  
    XX <- cbind(w, X)
    yhat <- predict(lm_fit, XX)
    return(as.vector(yhat))
  }
  
  return(list(
    model = lm_fit,
    tau_estimated = tau,
    generate_outcomes = ihdp_y_gen,
    design_matrix_formula = fml,
    lambda_used = best_lam
  ))
}

#' Run Single IHDP Analysis
#' 
#' Runs all randomization methods on the IHDP dataset once
#' 
#' @param X Covariate matrix
#' @param outcome_function Function to generate outcomes
#' @param true_tau True treatment effect
#' @param analysis_params Analysis parameters
#' 
#' @return Data frame with results
run_single_ihdp_analysis <- function(X, outcome_function, true_tau, analysis_params) {
  
  n <- nrow(X)
  methods_to_run <- c("Simple", "HCM", "HCM-ReR", "PWD-ReR")
  
  results_list <- list()
  runtimes <- list()
  
  # Simple Randomization
  start_time <- Sys.time()
  simple_result <- simple_randomization(n, seed = analysis_params$rep_seed)
  runtimes$Simple <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  results_list$Simple <- simple_result
  
  # Hierarchical Clustering Matching
  start_time <- Sys.time()
  hcm_result <- HCM(X, seed = analysis_params$rep_seed)
  runtimes$HCM <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  results_list$HCM <- hcm_result
  
  # HCM with Rerandomization
  start_time <- Sys.time()
  hcm_rer_result <- HCM_ReR(
    pa = analysis_params$pa, 
    X = X, 
    var_explained = analysis_params$var_explained,
    n_budget = analysis_params$n_budget, 
    seed = analysis_params$rep_seed
  )
  runtimes$`HCM-ReR` <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  results_list$`HCM-ReR` <- hcm_rer_result
  
  # PWD Rerandomization
  start_time <- Sys.time()
  pwd_rer_result <- PWD_ReR(
    pa = analysis_params$pa, 
    X = X, 
    var_explained = analysis_params$var_explained,
    n_budget = analysis_params$n_budget, 
    seed = analysis_params$rep_seed
  )
  runtimes$`PWD-ReR` <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  results_list$`PWD-ReR` <- pwd_rer_result
  
  # Evaluate each method
  evaluation_results <- data.frame()
  
  for (method_name in names(results_list)) {
    result <- results_list[[method_name]]
    runtime <- runtimes[[method_name]]
    
    # Generate outcomes
    y <- outcome_function(X, result$allocation)
    
    # Comprehensive evaluation
    evaluation <- comprehensive_evaluation(
      X = X, 
      result = result, 
      y = y, 
      tau_true = true_tau, 
      runtime_seconds = runtime
    )
    
    # Create evaluation row
    eval_row <- data.frame(
      replication = analysis_params$rep_id,
      method = method_name,
      sample_size = n,
      n_covariates = ncol(X),
      
      # Balance metrics
      mahalanobis_distance = evaluation$mahalanobis_distance,
      max_smd = evaluation$max_smd,
      mean_smd = evaluation$mean_smd,
      balance_score = evaluation$percent_smd_below_01,
      
      # Efficiency metrics
      iterations = evaluation$iterations_used,
      runtime_seconds = evaluation$runtime_seconds,
      budget_exhausted = evaluation$budget_exhausted,
      
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

#' Run IHDP Analysis Study
#' 
#' Runs multiple replications of the IHDP analysis
#' 
#' @param ihdp_data Preprocessed IHDP data
#' @param outcome_model Fitted outcome model
#' @param n_replications Number of replications
#' @param analysis_params Analysis parameters
#' @param save_results Whether to save results
#' @param results_dir Results directory
#' 
#' @return Analysis results
run_ihdp_analysis_study <- function(ihdp_data, outcome_model, n_replications = 1000,
                                    analysis_params = list(pa = 0.05, var_explained = 0.7, n_budget = 1000),
                                    save_results = TRUE, results_dir = "real_world_results") {
  
  cat("=== IHDP REAL-WORLD ANALYSIS STUDY ===\n")
  cat("Sample size:", nrow(ihdp_data$data), "\n")
  cat("Number of covariates:", ncol(ihdp_data$data) - 2, "\n")
  cat("Number of replications:", n_replications, "\n")
  cat("Estimated treatment effect:", round(outcome_model$tau_estimated, 4), "\n\n")
  
  # Create results directory
  if (save_results && !dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Prepare data
  ihdp <- ihdp_data$data
  n <- nrow(ihdp)
  
  # Create covariate matrix for analysis
  vars <- paste(paste('x', 1:(ncol(ihdp)-2), sep=''), collapse = '+')
  fml <- as.formula(paste('~', vars, '- 1'))
  X <- as.matrix(model.matrix(fml, data = ihdp))
  X <- scale(X)
  
  # Run replications
  all_results <- data.frame()
  start_time <- Sys.time()
  
  for (rep in 1:n_replications) {
    if (rep %% 100 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat("Completed", rep, "/", n_replications, "replications in", 
          round(elapsed, 2), "minutes\n")
    }
    
    # Set up replication parameters
    rep_params <- c(analysis_params, list(rep_id = rep, rep_seed = 2020 + rep))
    
    # Run single analysis
    tryCatch({
      rep_results <- run_single_ihdp_analysis(
        X = X,
        outcome_function = outcome_model$generate_outcomes,
        true_tau = outcome_model$tau_estimated,
        analysis_params = rep_params
      )
      
      all_results <- rbind(all_results, rep_results)
      
    }, error = function(e) {
      cat("Error in replication", rep, ":", e$message, "\n")
    })
  }
  
  # Calculate summary statistics
  summary_results <- all_results %>%
    group_by(method) %>%
    summarise(
      n_replications = n(),
      
      # Balance metrics
      mean_mahalanobis = mean(mahalanobis_distance, na.rm = TRUE),
      median_mahalanobis = median(mahalanobis_distance, na.rm = TRUE),
      sd_mahalanobis = sd(mahalanobis_distance, na.rm = TRUE),
      mean_max_smd = mean(max_smd, na.rm = TRUE),
      mean_balance_score = mean(balance_score, na.rm = TRUE),
      
      # Efficiency metrics
      mean_iterations = mean(iterations, na.rm = TRUE),
      median_iterations = median(iterations, na.rm = TRUE),
      mean_runtime = mean(runtime_seconds, na.rm = TRUE),
      budget_exhausted_rate = mean(budget_exhausted, na.rm = TRUE),
      
      # Treatment effect metrics
      mean_tau_hat = mean(tau_hat, na.rm = TRUE),
      sd_tau_hat = sd(tau_hat, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_abs_bias = mean(abs_bias, na.rm = TRUE),
      median_abs_bias = median(abs_bias, na.rm = TRUE),
      rmse = sqrt(mean(mse, na.rm = TRUE)),
      bias_sd = sd(bias, na.rm = TRUE),
      
      # Coverage probabilities (95% CI)
      ci_lower = quantile(tau_hat, 0.025, na.rm = TRUE),
      ci_upper = quantile(tau_hat, 0.975, na.rm = TRUE),
      
      .groups = 'drop'
    )
  
  # Add rankings
  summary_results$balance_rank <- rank(summary_results$mean_mahalanobis)
  summary_results$efficiency_rank <- rank(summary_results$mean_iterations)
  summary_results$bias_rank <- rank(summary_results$mean_abs_bias)
  summary_results$overall_rank <- rank(summary_results$mean_abs_bias + 
                                         0.5 * summary_results$balance_rank / nrow(summary_results) +
                                         0.3 * summary_results$efficiency_rank / nrow(summary_results))
  
  # Calculate runtime
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # Save results
  if (save_results) {
    write.csv(all_results, 
              file.path(results_dir, "ihdp_individual_results.csv"), 
              row.names = FALSE)
    write.csv(summary_results, 
              file.path(results_dir, "ihdp_summary_results.csv"), 
              row.names = FALSE)
    
    # Save analysis metadata
    metadata <- list(
      dataset = "IHDP",
      sample_size = nrow(ihdp_data$data),
      n_covariates = ncol(X),
      n_replications = n_replications,
      estimated_treatment_effect = outcome_model$tau_estimated,
      analysis_parameters = analysis_params,
      total_runtime_minutes = total_time
    )
    
    saveRDS(metadata, file.path(results_dir, "ihdp_analysis_metadata.rds"))
  }
  
  cat("IHDP analysis completed in", round(total_time, 2), "minutes\n\n")
  
  return(list(
    summary = summary_results,
    individual_results = all_results,
    metadata = list(
      sample_size = nrow(ihdp_data$data),
      n_covariates = ncol(X),
      estimated_treatment_effect = outcome_model$tau_estimated,
      total_runtime_minutes = total_time
    )
  ))
}

#' Generate IHDP Analysis Report
#' 
#' Generates comprehensive report for IHDP analysis results
#' 
#' @param analysis_results Results from run_ihdp_analysis_study
#' @param save_plots Whether to save plots
#' @param results_dir Results directory
#' 
#' @export
generate_ihdp_report <- function(analysis_results, save_plots = TRUE, 
                                 results_dir = "real_world_results") {
  
  if (!require(ggplot2)) stop("ggplot2 required for report generation")
  
  summary_results <- analysis_results$summary
  individual_results <- analysis_results$individual_results
  
  # Create plots directory
  plots_dir <- file.path(results_dir, "plots")
  if (save_plots && !dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Plot 1: Balance performance
  p1 <- ggplot(summary_results, aes(x = reorder(method, mean_mahalanobis), 
                                    y = mean_mahalanobis)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_mahalanobis - sd_mahalanobis/sqrt(n_replications),
                      ymax = mean_mahalanobis + sd_mahalanobis/sqrt(n_replications)),
                  width = 0.2) +
    scale_y_log10() +
    theme_minimal() +
    labs(title = "Balance Performance on IHDP Dataset",
         subtitle = "Mean Mahalanobis Distance (Lower is Better)",
         x = "Method", y = "Mahalanobis Distance (log scale)")
  
  if (save_plots) ggsave(file.path(plots_dir, "ihdp_balance_performance.png"), p1, width = 10, height = 6)
  print(p1)
  
  # Plot 2: Treatment effect estimation
  p2 <- ggplot(individual_results, aes(x = method, y = tau_hat)) +
    geom_boxplot(alpha = 0.7, fill = "lightcoral") +
    geom_hline(yintercept = analysis_results$metadata$estimated_treatment_effect, 
               linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Treatment Effect Estimation Distribution",
         subtitle = "Red line = True Effect",
         x = "Method", y = "Estimated Treatment Effect")
  
  if (save_plots) ggsave(file.path(plots_dir, "ihdp_treatment_effect_distribution.png"), p2, width = 10, height = 6)
  print(p2)
  
  # Plot 3: Bias comparison
  p3 <- ggplot(summary_results, aes(x = reorder(method, mean_abs_bias), y = mean_abs_bias)) +
    geom_col(fill = "lightgreen", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Treatment Effect Estimation Bias",
         subtitle = "Mean Absolute Bias (Lower is Better)",
         x = "Method", y = "Mean Absolute Bias")
  
  if (save_plots) ggsave(file.path(plots_dir, "ihdp_bias_comparison.png"), p3, width = 10, height = 6)
  print(p3)
  
  # Plot 4: Efficiency vs Accuracy trade-off
  p4 <- ggplot(summary_results, aes(x = mean_iterations, y = mean_abs_bias, 
                                    color = method, size = mean_mahalanobis)) +
    geom_point(alpha = 0.8) +
    geom_text(aes(label = method), vjust = -0.5, hjust = 0.5, size = 3) +
    scale_x_log10() +
    theme_minimal() +
    labs(title = "Efficiency vs Accuracy Trade-off (IHDP Dataset)",
         subtitle = "Lower-left is optimal (low iterations, low bias)",
         x = "Mean Iterations (log scale)", y = "Mean Absolute Bias",
         color = "Method", size = "Balance (Mahalanobis)")
  
  if (save_plots) ggsave(file.path(plots_dir, "ihdp_efficiency_accuracy_tradeoff.png"), p4, width = 12, height = 8)
  print(p4)
  
  # Generate text report
  cat("\n=== IHDP REAL-WORLD ANALYSIS REPORT ===\n\n")
  
  cat("Dataset Information:\n")
  cat("- Sample size:", analysis_results$metadata$sample_size, "\n")
  cat("- Number of covariates:", analysis_results$metadata$n_covariates, "\n")
  cat("- Estimated treatment effect:", round(analysis_results$metadata$estimated_treatment_effect, 4), "\n\n")
  
  cat("Method Performance Rankings:\n")
  rankings <- summary_results[order(summary_results$overall_rank), 
                              c("method", "overall_rank", "balance_rank", "efficiency_rank", "bias_rank")]
  print(rankings)
  
  cat("\n=== KEY FINDINGS ===\n")
  best_overall <- summary_results$method[which.min(summary_results$overall_rank)]
  best_balance <- summary_results$method[which.min(summary_results$mean_mahalanobis)]
  best_efficiency <- summary_results$method[which.min(summary_results$mean_iterations)]
  best_bias <- summary_results$method[which.min(summary_results$mean_abs_bias)]
  
  cat("Best Overall Method:", best_overall, "\n")
  cat("Best Balance Method:", best_balance, "\n")
  cat("Most Efficient Method:", best_efficiency, "\n")
  cat("Least Biased Method:", best_bias, "\n\n")
  
  # Method-specific insights
  if (best_overall == "HCM-ReR") {
    cat("HCM-ReR demonstrates superior performance by combining hierarchical clustering\n")
    cat("with PCA-weighted rerandomization, achieving excellent balance while maintaining\n")
    cat("computational efficiency.\n")
  } else if (best_overall == "PWD-ReR") {
    cat("PWD-ReR shows optimal performance through its focus on the most important\n")
    cat("covariate dimensions identified by PCA, leading to efficient balance optimization.\n")
  } else if (best_overall == "HCM") {
    cat("HCM provides the best performance through its simple but effective approach\n")
    cat("of matching similar units via hierarchical clustering.\n")
  }
  
  cat("\n=== STATISTICAL SUMMARY ===\n")
  detailed_summary <- summary_results[, c("method", "mean_mahalanobis", "mean_iterations", 
                                          "mean_abs_bias", "rmse", "budget_exhausted_rate")]
  print(round(detailed_summary, 4))
  
  invisible(list(
    plots = list(balance = p1, treatment_effects = p2, bias = p3, tradeoff = p4),
    best_methods = list(overall = best_overall, balance = best_balance, 
                        efficiency = best_efficiency, bias = best_bias)
  ))
}

#' Main IHDP Analysis Function
#' 
#' Complete workflow for IHDP analysis
#' 
#' @param data_path Path to IHDP data files
#' @param n_replications Number of analysis replications
#' @param analysis_params Analysis parameters
#' @param save_results Whether to save results
#' @param results_dir Results directory
#' 
#' @return Complete analysis results
main_ihdp_analysis <- function(data_path = "data/", n_replications = 1000,
                               analysis_params = list(pa = 0.05, var_explained = 0.7, n_budget = 1000),
                               save_results = TRUE, results_dir = "real_world_results") {
  
  cat("=== STARTING IHDP REAL-WORLD ANALYSIS ===\n\n")
  
  # Load and preprocess data
  ihdp_data <- load_and_preprocess_ihdp(data_path)
  
  # Fit outcome model
  outcome_model <- fit_ihdp_outcome_model(ihdp_data, use_interactions = TRUE)
  
  # Run analysis study
  analysis_results <- run_ihdp_analysis_study(
    ihdp_data = ihdp_data,
    outcome_model = outcome_model,
    n_replications = n_replications,
    analysis_params = analysis_params,
    save_results = save_results,
    results_dir = results_dir
  )
  
  # Generate comprehensive report
  report <- generate_ihdp_report(
    analysis_results = analysis_results,
    save_plots = save_results,
    results_dir = results_dir
  )
  
  cat("\n=== IHDP ANALYSIS COMPLETED ===\n")
  if (save_results) {
    cat("Results saved to:", results_dir, "\n")
  }
  
  return(list(
    data = ihdp_data,
    outcome_model = outcome_model,
    analysis_results = analysis_results,
    report = report
  ))
}

# Uncomment to run the complete IHDP analysis
# ihdp_results <- main_ihdp_analysis()