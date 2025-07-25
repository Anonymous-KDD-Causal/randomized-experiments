# Randomized Experiments: HCM, HCM-ReR, and PWD-ReR Methods

[![R-CMD-check](https://github.com/yourusername/randomized-experiments/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/randomized-experiments/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This R package implements three advanced methods for generating balanced treatment allocations in randomized experiments:

- **HCM**: Hierarchical Clustering Matched - Basic pair-based randomization using hierarchical clustering
- **HCM-ReR**: Hierarchical Clustering Matched Re-Randomization - Enhanced HCM with PCA weighting and re-randomization
- **PWD-ReR**: Principal Component Weighted Distance Re-Randomization - Complete re-randomization with PCA-weighted balance criteria

## Installation

```r
# Install from GitHub
devtools::install_github("yourusername/randomized-experiments-main")

library(randomizedexperiments)

# Generate example data
X <- generate_covariates(n = 100, p = 10, rho = 0.3, seed = 123)

# Run all three methods
result_hcm <- hcm(X, seed = 123)
result_hcm_rer <- hcm_rer(X, pa = 0.1, max_iterations = 500, seed = 123)
result_pwd_rer <- pwd_rer(X, pa = 0.1, max_iterations = 500, seed = 123)

# Compare methods using built-in function
comparison <- compare_methods(X, seed = 123)
print(comparison$summary)

# Evaluate balance
balance_hcm <- calculate_balance_metrics(X, result_hcm$allocation)
balance_hcm_rer <- calculate_balance_metrics(X, result_hcm_rer$allocation)
balance_pwd_rer <- calculate_balance_metrics(X, result_pwd_rer$allocation)

print(paste("HCM Balance Score:", round(balance_hcm$balance_score, 3)))
print(paste("HCM-ReR Balance Score:", round(balance_hcm_rer$balance_score, 3)))
print(paste("PWD-ReR Balance Score:", round(balance_pwd_rer$balance_score, 3)))
