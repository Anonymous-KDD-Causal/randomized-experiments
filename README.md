# Hierarchical Clustering Matching and Rerandomization for Causality

[![R-CMD-check](https://github.com/yourusername/randomized-experiments/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/randomized-experiments/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This R package implements three advanced methods for generating balanced treatment allocations in randomized experiments:

- **HCM**: Hierarchical Clustering Matching - Basic pair-based randomization using hierarchical clustering
- **HCM-ReR**: Hierarchical Clustering Matching Re-Randomization - An enhanced version of HCM. It uses hierarchical clustering to form pairs and re-randomizes treatment assignments within pairs 
- **PWD-ReR**: Principal Component Weighted Distance Re-Randomization - Complete re-randomization with PCA-weighted balance criteria

## Installation

```r
# Install from GitHub
devtools::install_github("yourusername/randomized-experiments-main")

# Quick Start

source("R/hcm.R")
source("R/hcm_rer.R")
source("R/pwd_rer.R")
source("R/data_generation.R")

seed = 2000

# Generate synthetic experimental data
X <- generate_covariates(n = 200, d = 50, rho = 0.2, seed = seed)
X <- scale(X)  # Always standardize

# HCM
hcm_result <- hcm(X, seed = seed)

# HCM-ReR
hcm_rer_result <- hcm_rer(X, pa = 0.05, n_budget = 500, seed = seed)

# PWD-ReR
pwd_rer_result <- pwd_rer(X, pa = 0.05, var_explained = 0.7, seed = seed)

# Evaluate all methods
eval_hcm <- evaluate_hcm(X, hcm_result)
eval_hcm_rer <- evaluate_hcm_rer(X, hcm_rer_result)  
eval_pwd_rer <- evaluate_pwd_rer(X, pwd_rer_result)

# ============================================================================
# EXTRACT RESULTS AND CALCULATE BIAS
# ============================================================================

# Set parameters for treatment effect estimation
tau_true <- 1.0
p <- ncol(X)
beta <- rep(0.5, p)
seed_bias <- seed

# Calculate treatment effects for each method
# HCM
set.seed(seed_bias)
y_hcm <- generate_response(X, hcm_result$allocation, beta, tau = tau_true, seed = seed_bias)
tau_hat_hcm <- mean(y_hcm[hcm_result$allocation == 1]) - mean(y_hcm[hcm_result$allocation == 0])

# HCM-ReR
set.seed(seed_bias + 1)
y_hcm_rer <- generate_response(X, hcm_rer_result$allocation, beta, tau = tau_true, seed = seed_bias + 1)
tau_hat_hcm_rer <- mean(y_hcm_rer[hcm_rer_result$allocation == 1]) - mean(y_hcm_rer[hcm_rer_result$allocation == 0])

# PWD-ReR
set.seed(seed_bias + 2)
y_pwd_rer <- generate_response(X, pwd_rer_result$allocation, beta, tau = tau_true, seed = seed_bias + 2)
tau_hat_pwd_rer <- mean(y_pwd_rer[pwd_rer_result$allocation == 1]) - mean(y_pwd_rer[pwd_rer_result$allocation == 0])

# ============================================================================
# CREATE RESULTS TABLE
# ============================================================================

results_table <- data.frame(
  Method = c("HCM", "HCM-ReR", "PWD-ReR"),
  Mahalanobis_Distance = c(eval_hcm$mahalanobis_distance, 
                           eval_hcm_rer$mahalanobis_distance, 
                           eval_pwd_rer$mahalanobis_distance),
  Max_SMD = c(eval_hcm$max_smd, 
              eval_hcm_rer$max_smd, 
              eval_pwd_rer$max_smd),
  Mean_SMD = c(eval_hcm$mean_smd, 
               eval_hcm_rer$mean_smd, 
               eval_pwd_rer$mean_smd),
  Balance_Score = c(eval_hcm$balance_score, 
                    eval_hcm_rer$balance_score, 
                    eval_pwd_rer$balance_score),
  True_Tau = rep(tau_true, 3),
  Estimated_Tau = c(tau_hat_hcm, tau_hat_hcm_rer, tau_hat_pwd_rer),
  Bias = c(tau_hat_hcm - tau_true, 
           tau_hat_hcm_rer - tau_true, 
           tau_hat_pwd_rer - tau_true),
  Abs_Bias = c(abs(tau_hat_hcm - tau_true), 
               abs(tau_hat_hcm_rer - tau_true), 
               abs(tau_hat_pwd_rer - tau_true))
)

# Round numeric columns
numeric_cols <- sapply(results_table, is.numeric)
results_table[numeric_cols] <- round(results_table[numeric_cols], 4)

# Print results table
cat("=== METHOD RESULTS TABLE ===\n")
print(results_table, row.names = FALSE)


