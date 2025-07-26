# Hierarchical Clustering Matching and Rerandomization for Causality

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
