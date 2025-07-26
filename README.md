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

library(devtools)

# Install from local directory
devtools::install_local("yourusername/randomized-experiments-main",force=TRUE)

library(randomizedexperiments)
seed = 2025

X <- generate_covariates(n = 500, p = 50, rho = 0.2, seed = seed)

# Run all three methods
result_hcm <- hcm(X, seed = seed)
result_hcm_rer <- hcm_rer(X, pa = 0.05, max_iterations = 1000, seed = seed)
result_pwd_rer <- pwd_rer(X, pa = 0.05, max_iterations = 1000, seed = seed)

# Compare methods
comparison <- compare_methods(X, seed = seed, pa = 0.05, 
                              true_tau = 1,sigma = 0.05, 
                              response_type = "linear")
print(comparison$summary)

