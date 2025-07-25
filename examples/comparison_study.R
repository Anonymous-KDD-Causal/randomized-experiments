#' Comprehensive Comparison Study: HCM vs HCM-ReR vs PWD-ReR
#' 
#' This example demonstrates all three methods and provides detailed comparison.

library(randomizedexperiments)

# Set parameters
n <- 200          # Sample size (must be even for HCM methods)
p <- 15           # Number of covariates  
rho <- 0.4        # Correlation between covariates
pa <- 0.15        # Acceptance probability for re-randomization
max_iter <- 1000  # Maximum iterations
seed <- 54321     # Random seed

cat("=== Comprehensive Methods Comparison ===\n")
cat("Methods: HCM, HCM-ReR, PWD-ReR\n")
cat("Parameters:\n")
cat("- Sample size (n):", n, "\n")
cat("- Covariates (p):", p, "\n") 
cat("- Correlation (rho):", rho, "\n")
cat("- Acceptance prob (pa):", pa, "\n\n")

# Generate covariates
cat("Generating covariates...\n")
X <- generate_covariates(n = n, p = p, rho = rho, seed = seed)

# Method 1: Basic HCM
cat("\n=== Running HCM (Basic Hierarchical Clustering Matched) ===\n")
start_time <- Sys.time()
result_hcm <- hcm(X, seed = seed)
time_hcm <- as.numeric(Sys.time() - start_time)

# Method 2: HCM-ReR
cat("\n=== Running HCM-ReR (with Re-Randomization) ===\n")
start_time <- Sys.time()
result_hcm_rer <- hcm_rer(X, pa = pa, max_iterations = max_iter, seed = seed)
time_hcm_rer <- as.numeric(Sys.time() - start_time)

# Method 3: PWD-ReR  
cat("\n=== Running PWD-ReR (Complete Re-Randomization) ===\n")
start_time <- Sys.time()
result_pwd_rer <- pwd_rer(X, pa = pa, max_iterations = max_iter, seed = seed)
time_pwd_rer <- as.numeric(Sys.time() - start_time)

# Evaluate all methods
cat("\n=== Evaluating Balance Metrics ===\n")
eval_hcm <- evaluate_allocation(X, result_hcm$allocation, true_tau = 1, seed = seed)
eval_hcm_rer <- evaluate_allocation(X, result_hcm_rer$allocation, true_tau = 1, seed = seed)
eval_pwd_rer <- evaluate_allocation(X, result_pwd_rer$allocation, true_tau = 1, seed = seed)

# Print detailed comparison
cat("\n=== DETAILED RESULTS COMPARISON ===\n")
cat(sprintf("%-25s %-10s %-15s %-10s\n", "Metric", "HCM", "HCM-ReR", "PWD-ReR"))
cat(sprintf("%-25s %-10s %-15s %-10s\n", 
            "-------------------------", "----------", "---------------", "----------"))

# Basic information
cat(sprintf("%-25s %-10s %-15s %-10s\n", "Method Type", "Matching", "Match+ReR", "Complete ReR"))
cat(sprintf("%-25s %-10s %-15s %-10s\n", "Requires Even N", "Yes", "Yes", "No"))

# Acceptance and iterations
cat(sprintf("%-25s %-10s %-15s %-10s\n", 
            "Accepted", "N/A", 
            result_hcm_rer$accepted, 
            result_pwd_rer$accepted))
cat(sprintf("%-25s %-10d %-15d %-10d\n", 
            "Iterations", 1, 
            result_hcm_rer$iterations, 
            result_pwd_rer$iterations))

# Runtime
cat(sprintf("%-25s %-10.2f %-15.2f %-10.2f\n", 
            "Runtime (seconds)", time_hcm, time_hcm_rer, time_pwd_rer))

# Balance metrics
cat(sprintf("%-25s %-10.3f %-15.3f %-10.3f\n", 
            "Balance Score", eval_hcm$balance_score, 
            eval_hcm_rer$balance_score, eval_pwd_rer$balance_score))
cat(sprintf("%-25s %-10.3f %-15.3f %-10.3f\n", 
            "Max SMD", eval_hcm$max_smd, 
            eval_hcm_rer$max_smd, eval_pwd_rer$max_smd))
cat(sprintf("%-25s %-10.3f %-15.3f %-10.3f\n", 
            "Mean SMD", eval_hcm$mean_smd, 
            eval_hcm_rer$mean_smd, eval_pwd_rer$mean_smd))
cat(sprintf("%-25s %-10.1f %-15.1f %-10.1f\n", 
            "Mahalanobis Distance", eval_hcm$mahalanobis_distance, 
            eval_hcm_rer$mahalanobis_distance, eval_pwd_rer$mahalanobis_distance))

# Treatment effect estimation
cat(sprintf("%-25s %-10.4f %-15.4f %-10.4f\n", 
            "Treatment Effect Est", eval_hcm$tau_hat, 
            eval_hcm_rer$tau_hat, eval_pwd_rer$tau_hat))
cat(sprintf("%-25s %-10.4f %-15.4f %-10.4f\n", 
            "Absolute Bias", eval_hcm$abs_bias, 
            eval_hcm_rer$abs_bias, eval_pwd_rer$abs_bias))

# Sample size efficiency
cat(sprintf("%-25s %-10.3f %-15.3f %-10.3f\n", 
            "Sample Size Efficiency", eval_hcm$sample_size_efficiency, 
            eval_hcm_rer$sample_size_efficiency, eval_pwd_rer$sample_size_efficiency))

# Analysis of matched pairs (HCM methods only)
cat("\n=== PAIR ANALYSIS (HCM Methods) ===\n")

# HCM pair analysis
hcm_pairs <- analyze_pair_distances(result_hcm, X)
cat("HCM Pairs:\n")
cat("  Number of pairs:", hcm_pairs$n_pairs, "\n")
cat("  Mean within-pair distance:", sprintf("%.3f", hcm_pairs$mean_distance), "\n")
cat("  SD of within-pair distances:", sprintf("%.3f", hcm_pairs$sd_distance), "\n")

# HCM-ReR pair analysis
hcm_rer_pairs <- analyze_pair_distances(result_hcm_rer, X)
cat("\nHCM-ReR Pairs:\n")
cat("  Number of pairs:", hcm_rer_pairs$n_pairs, "\n")
cat("  Mean within-pair distance:", sprintf("%.3f", hcm_rer_pairs$mean_distance), "\n")
cat("  SD of within-pair distances:", sprintf("%.3f", hcm_rer_pairs$sd_distance), "\n")

# PCA analysis for HCM-ReR and PWD-ReR
cat("\n=== PCA ANALYSIS (Re-Randomization Methods) ===\n")
cat("HCM-ReR:\n")
cat("  Components used:", result_hcm_rer$k_components, "\n")
cat("  Top 3 variable weights:", sprintf("%.3f", head(sort(result_hcm_rer$var_weights, decreasing = TRUE), 3)), "\n")

cat("\nPWD-ReR:\n")
cat("  Components used:", result_pwd_rer$k_components, "\n")
cat("  Variance explained:", sprintf("%.3f", result_pwd_rer$pca_variance_explained), "\n")
cat("  Top 3 variable weights:", sprintf("%.3f", head(sort(result_pwd_rer$var_weights, decreasing = TRUE), 3)), "\n")

# Use compare_methods function for additional analysis
cat("\n=== USING BUILT-IN COMPARISON FUNCTION ===\n")
comparison <- compare_methods(X, pa = pa, max_iterations = max_iter, seed = seed)
print(comparison$summary)

cat("\nBest performing method:", comparison$best_method, "\n")

# Optional: Create plots if ggplot2 is available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("\n=== Creating Comparison Plots ===\n")
  
  # Plot method comparison
  p1 <- plot_method_comparison(comparison)
  print(p1)
  
  # Balance score evolution for re-randomization methods
  library(ggplot2)
  
  # Create balance history plot
  if (length(result_hcm_rer$balance_history) > 0 && length(result_pwd_rer$balance_history) > 0) {
    balance_data <- data.frame(
      Iteration = c(1:length(result_hcm_rer$balance_history), 
                    1:length(result_pwd_rer$balance_history)),
      Balance_Score = c(result_hcm_rer$balance_history, result_pwd_rer$balance_history),
      Method = c(rep("HCM-ReR", length(result_hcm_rer$balance_history)),
                 rep("PWD-ReR", length(result_pwd_rer$balance_history)))
    )
    
    p2 <- ggplot(balance_data, aes(x = Iteration, y = Balance_Score, color = Method)) +
      geom_line(alpha = 0.7) +
      geom_hline(yintercept = result_hcm_rer$balance_threshold, linetype = "dashed", color = "red", alpha = 0.7) +
      theme_minimal() +
      labs(title = "Balance Score Evolution During Re-Randomization",
           x = "Iteration",
           y = "Weighted Balance Score",
           subtitle = "Dashed line shows acceptance threshold") +
      scale_color_brewer(type = "qual", palette = "Set1")
    
    print(p2)
  }
}

cat("\n=== RECOMMENDATIONS ===\n")
cat("Based on the results:\n")

best_balance <- which.max(c(eval_hcm$balance_score, eval_hcm_rer$balance_score, eval_pwd_rer$balance_score))
methods_names <- c("HCM", "HCM-ReR", "PWD-ReR")

cat("- Best balance achieved by:", methods_names[best_balance], "\n")

if (result_hcm_rer$accepted && result_pwd_rer$accepted) {
  cat("- Both re-randomization methods found acceptable allocations\n")
} else if (result_hcm_rer$accepted) {
  cat("- Only HCM-ReR found acceptable allocation\n")
} else if (result_pwd_rer$accepted) {
  cat("- Only PWD-ReR found acceptable allocation\n")
} else {
  cat("- Neither re-randomization method found acceptable allocation within iteration limit\n")
  cat("- Consider increasing max_iterations or pa threshold\n")
}

if (time_hcm_rer > 2 * time_hcm && eval_hcm_rer$balance_score - eval_hcm$balance_score < 0.1) {
  cat("- HCM-ReR provides marginal improvement over HCM at significant computational cost\n")
}

cat("\nComparison study completed successfully!\n")