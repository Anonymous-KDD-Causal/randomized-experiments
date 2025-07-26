#' Data Generation Functions for Randomization Method Testing
#' 
#' This file contains functions for generating synthetic datasets with various
#' correlation structures and outcome models for testing randomization methods.

#' @import MASS
#' @import stats

#' Generate Correlated Covariates
#' 
#' Generates a matrix of correlated covariates using multivariate normal distribution.
#' 
#' @param n Sample size
#' @param p Number of covariates
#' @param correlation_structure Type of correlation structure: "exchangeable", "ar1", "block", or "custom"
#' @param rho Correlation parameter(s)
#' @param seed Random seed
#' 
#' @return Matrix of covariates (n x p)
#' @export
generate_covariates <- function(n, p, correlation_structure = "exchangeable", 
                                rho = 0.5, seed = 2020) {
  set.seed(seed)
  
  # Create correlation matrix based on structure
  if (correlation_structure == "exchangeable") {
    # All pairwise correlations equal to rho
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
    
  } else if (correlation_structure == "ar1") {
    # AR(1) structure: Sigma[i,j] = rho^|i-j|
    Sigma <- rho^abs(outer(1:p, 1:p, "-"))
    
  } else if (correlation_structure == "block") {
    # Block diagonal structure
    if (length(rho) == 1) rho <- rep(rho, ceiling(p/5))
    block_size <- ceiling(p / length(rho))
    Sigma <- diag(p)
    
    for (i in seq_along(rho)) {
      start_idx <- (i-1) * block_size + 1
      end_idx <- min(i * block_size, p)
      if (start_idx <= p) {
        block_indices <- start_idx:end_idx
        Sigma[block_indices, block_indices] <- rho[i]
        diag(Sigma[block_indices, block_indices]) <- 1
      }
    }
    
  } else if (correlation_structure == "custom") {
    # Require rho to be a full correlation matrix
    if (!is.matrix(rho) || nrow(rho) != p || ncol(rho) != p) {
      stop("For custom correlation structure, rho must be a p x p correlation matrix")
    }
    Sigma <- rho
    
  } else {
    stop("Unknown correlation structure. Use 'exchangeable', 'ar1', 'block', or 'custom'")
  }
  
  # Generate multivariate normal data
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  return(X)
}

#' Generate Factor Structure Data
#' 
#' Generates covariates with underlying factor structure.
#' 
#' @param n Sample size
#' @param n_factors Number of underlying factors
#' @param vars_per_factor Number of variables per factor
#' @param noise_vars Number of pure noise variables
#' @param factor_loading Loading of variables on factors
#' @param noise_level Level of noise in the factor model
#' @param seed Random seed
#' 
#' @return Matrix of covariates (n x p)
#' @export
generate_factor_data <- function(n, n_factors = 3, vars_per_factor = 5, noise_vars = 5,
                                 factor_loading = 0.8, noise_level = 0.3, seed = 2020) {
  set.seed(seed)
  
  # Generate factors
  factors <- matrix(rnorm(n * n_factors), n, n_factors)
  
  # Generate factor-related variables
  X_factors <- matrix(0, n, n_factors * vars_per_factor)
  
  for (f in 1:n_factors) {
    start_col <- (f-1) * vars_per_factor + 1
    end_col <- f * vars_per_factor
    
    for (v in start_col:end_col) {
      X_factors[, v] <- factor_loading * factors[, f] + 
        sqrt(1 - factor_loading^2) * rnorm(n)
    }
  }
  
  # Generate noise variables
  if (noise_vars > 0) {
    X_noise <- matrix(rnorm(n * noise_vars), n, noise_vars)
    X <- cbind(X_factors, X_noise)
  } else {
    X <- X_factors
  }
  
  # Add column names
  factor_names <- paste0("F", rep(1:n_factors, each = vars_per_factor), 
                         "V", rep(1:vars_per_factor, n_factors))
  if (noise_vars > 0) {
    noise_names <- paste0("N", 1:noise_vars)
    colnames(X) <- c(factor_names, noise_names)
  } else {
    colnames(X) <- factor_names
  }
  
  return(X)
}

#' Generate Treatment Effects with Heterogeneity
#' 
#' Generates outcomes with treatment effect heterogeneity based on covariates.
#' 
#' @param X Covariate matrix
#' @param w Treatment allocation vector
#' @param effect_type Type of treatment effect: "constant", "linear", "nonlinear", "interaction"
#' @param base_effect Base treatment effect size
#' @param effect_coefficients Coefficients for effect heterogeneity (if NULL, will be generated)
#' @param noise_sd Standard deviation of noise term
#' @param seed Random seed
#' 
#' @return Vector of outcomes
#' @export
generate_outcomes <- function(X, w, effect_type = "linear", base_effect = 1.0,
                              effect_coefficients = NULL, noise_sd = 1.0, seed = 2020) {
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  
  # Generate baseline outcome (control potential outcome)
  if (is.null(effect_coefficients)) {
    # Default: random coefficients with decreasing importance
    effect_coefficients <- rnorm(p) * exp(-seq(0, 2, length.out = p))
  }
  
  # Baseline outcome (what would happen without treatment)
  y0 <- X %*% effect_coefficients + rnorm(n, sd = noise_sd)
  
  # Generate treatment effects
  if (effect_type == "constant") {
    # Constant treatment effect
    treatment_effect <- rep(base_effect, n)
    
  } else if (effect_type == "linear") {
    # Linear treatment effect heterogeneity
    te_coefficients <- effect_coefficients * 0.3  # Smaller coefficients for treatment effect
    treatment_effect <- base_effect + X %*% te_coefficients
    
  } else if (effect_type == "nonlinear") {
    # Nonlinear treatment effect
    linear_combo <- X %*% effect_coefficients[1:min(3, p)]
    treatment_effect <- base_effect + 0.5 * tanh(linear_combo / p)
    
  } else if (effect_type == "interaction") {
    # Treatment effect depends on interactions
    if (p >= 2) {
      interaction_term <- X[, 1] * X[, 2]
      treatment_effect <- base_effect + 0.3 * interaction_term
    } else {
      treatment_effect <- rep(base_effect, n)
    }
    
  } else {
    stop("Unknown effect_type. Use 'constant', 'linear', 'nonlinear', or 'interaction'")
  }
  
  # Generate observed outcomes
  y1 <- y0 + treatment_effect  # Treated potential outcome
  y_observed <- w * y1 + (1 - w) * y0
  
  # Return results with additional information
  return(list(
    y_observed = y_observed,
    y0 = y0,
    y1 = y1,
    treatment_effect = treatment_effect,
    average_treatment_effect = mean(treatment_effect),
    effect_coefficients = effect_coefficients
  ))
}

#' Generate Complete Simulation Dataset
#' 
#' Generates a complete dataset for randomization method testing.
#' 
#' @param n Sample size
#' @param p Number of covariates
#' @param correlation_structure Correlation structure for covariates
#' @param rho Correlation parameter
#' @param effect_type Type of treatment effect
#' @param base_effect Base treatment effect
#' @param noise_sd Outcome noise standard deviation
#' @param seed Random seed
#' 
#' @return List containing covariates, true treatment effects, and data generation parameters
#' @export
generate_simulation_data <- function(n = 200, p = 20, correlation_structure = "ar1", 
                                     rho = 0.5, effect_type = "linear", 
                                     base_effect = 1.0, noise_sd = 1.0, seed = 2020) {
  set.seed(seed)
  
  # Generate covariates
  X <- generate_covariates(n, p, correlation_structure, rho, seed)
  
  # Standardize covariates
  X <- scale(X)
  
  # Generate a reference treatment allocation for outcome generation
  w_ref <- sample(c(rep(0, n/2), rep(1, n/2)))
  
  # Generate outcomes
  outcome_data <- generate_outcomes(X, w_ref, effect_type, base_effect, 
                                    effect_coefficients = NULL, noise_sd, seed)
  
  return(list(
    X = X,
    true_ate = outcome_data$average_treatment_effect,
    individual_effects = outcome_data$treatment_effect,
    effect_coefficients = outcome_data$effect_coefficients,
    parameters = list(
      n = n,
      p = p,
      correlation_structure = correlation_structure,
      rho = rho,
      effect_type = effect_type,
      base_effect = base_effect,
      noise_sd = noise_sd,
      seed = seed
    ),
    generate_outcomes_function = function(allocation) {
      generate_outcomes(X, allocation, effect_type, base_effect, 
                        outcome_data$effect_coefficients, noise_sd, seed + 1000)$y_observed
    }
  ))
}

#' Generate IHDP-like Dataset
#' 
#' Generates a dataset similar to the IHDP (Infant Health and Development Program) dataset.
#' 
#' @param n Sample size
#' @param seed Random seed
#' 
#' @return List containing covariates and outcome generation function
#' @export
generate_ihdp_like_data <- function(n = 747, seed = 2020) {
  set.seed(seed)
  
  # Continuous covariates (6 variables similar to IHDP)
  birth_weight <- rnorm(n, mean = 2500, sd = 500)  # Birth weight in grams
  head_circumference <- rnorm(n, mean = 33, sd = 2)  # Head circumference in cm
  gestational_age <- pmax(24, rnorm(n, mean = 38, sd = 3))  # Gestational age in weeks
  birth_order <- rpois(n, lambda = 1.5) + 1  # Birth order
  neonatal_health <- rnorm(n, mean = 0, sd = 1)  # Neonatal health index
  mother_age <- pmax(15, rnorm(n, mean = 25, sd = 6))  # Mother's age
  
  # Binary covariates (similar to IHDP categorical variables)
  sex <- rbinom(n, 1, 0.5)  # Child sex (0 = female, 1 = male)
  twin <- rbinom(n, 1, 0.05)  # Twin birth
  mother_married <- rbinom(n, 1, 0.6)  # Mother married
  mother_education_low <- rbinom(n, 1, 0.3)  # Mother education < high school
  mother_education_high <- rbinom(n, 1, 0.4)  # Mother education = high school
  smoking <- rbinom(n, 1, 0.2)  # Smoking during pregnancy
  firstborn <- rbinom(n, 1, 0.35)  # First born child
  alcohol <- rbinom(n, 1, 0.15)  # Alcohol during pregnancy
  drugs <- rbinom(n, 1, 0.05)  # Drug use during pregnancy
  work_during_pregnancy <- rbinom(n, 1, 0.7)  # Worked during pregnancy
  prenatal_care <- rbinom(n, 1, 0.9)  # Received prenatal care
  
  # Site indicators (11 sites in original IHDP)
  site <- sample(1:8, n, replace = TRUE, prob = rep(1/8, 8))
  site_matrix <- model.matrix(~ factor(site) - 1)
  colnames(site_matrix) <- paste0("site_", 1:8)
  
  # Combine all covariates
  X_continuous <- cbind(birth_weight, head_circumference, gestational_age, 
                        birth_order, neonatal_health, mother_age)
  X_binary <- cbind(sex, twin, mother_married, mother_education_low, mother_education_high,
                    smoking, firstborn, alcohol, drugs, work_during_pregnancy, prenatal_care)
  
  X <- cbind(X_continuous, X_binary, site_matrix)
  
  # Standardize continuous variables
  X[, 1:6] <- scale(X[, 1:6])
  
  # Create realistic effect structure
  # Treatment effect should depend on baseline characteristics
  effect_coefficients <- c(
    # Continuous variables
    0.3, -0.2, 0.1, -0.1, 0.4, -0.1,
    # Binary variables  
    0.1, -0.3, 0.2, -0.4, 0.3, -0.3, 0.2, -0.2, -0.4, 0.1, 0.2,
    # Site effects
    rep(0.1, 8)
  )
  
  # Outcome generation function
  generate_outcomes_function <- function(allocation) {
    # Baseline outcome
    y0 <- X %*% effect_coefficients + rnorm(n, sd = 1.2)
    
    # Treatment effect (positive but heterogeneous)
    treatment_effect <- 4.0 + 0.3 * X[, 1] - 0.2 * X[, 4] + 0.4 * X[, 5]  # Depends on birth weight, birth order, neonatal health
    
    # Observed outcome
    y_observed <- y0 + allocation * treatment_effect + rnorm(n, sd = 0.5)
    
    return(y_observed)
  }
  
  colnames(X) <- c("birth_weight", "head_circumference", "gestational_age", 
                   "birth_order", "neonatal_health", "mother_age",
                   "sex", "twin", "mother_married", "mother_education_low", "mother_education_high",
                   "smoking", "firstborn", "alcohol", "drugs", "work_during_pregnancy", "prenatal_care",
                   paste0("site_", 1:8))
  
  return(list(
    X = X,
    true_ate = 4.0,  # Average treatment effect
    generate_outcomes_function = generate_outcomes_function,
    variable_descriptions = list(
      continuous = c("birth_weight", "head_circumference", "gestational_age", "birth_order", "neonatal_health", "mother_age"),
      binary = c("sex", "twin", "mother_married", "mother_education_low", "mother_education_high", "smoking", "firstborn", "alcohol", "drugs", "work_during_pregnancy", "prenatal_care"),
      categorical = paste0("site_", 1:8)
    )
  ))
}

#' Generate Multiple Datasets for Simulation Study
#' 
#' Generates multiple datasets with different characteristics for comprehensive testing.
#' 
#' @param n_datasets Number of datasets to generate
#' @param ... Additional parameters passed to generate_simulation_data
#' 
#' @return List of generated datasets
#' @export
generate_multiple_datasets <- function(n_datasets = 5, ...) {
  datasets <- list()
  
  # Different configurations
  configurations <- list(
    list(correlation_structure = "exchangeable", rho = 0.3, p = 10),
    list(correlation_structure = "ar1", rho = 0.5, p = 15),
    list(correlation_structure = "ar1", rho = 0.8, p = 20),
    list(correlation_structure = "block", rho = c(0.3, 0.7, 0.5), p = 25),
    list(correlation_structure = "exchangeable", rho = 0.6, p = 30)
  )
  
  for (i in 1:min(n_datasets, length(configurations))) {
    config <- configurations[[i]]
    datasets[[i]] <- do.call(generate_simulation_data, c(config, list(...)))
    names(datasets)[i] <- paste0("dataset_", i, "_", config$correlation_structure, "_rho", config$rho)
  }
  
  return(datasets)
}

#' Summary of Generated Data
#' 
#' Provides a summary of generated dataset characteristics.
#' 
#' @param data_object Dataset object from generation functions
#' 
#' @export
summarize_generated_data <- function(data_object) {
  cat("=== Generated Data Summary ===\n")
  cat("Sample size (n):", nrow(data_object$X), "\n")
  cat("Number of covariates (p):", ncol(data_object$X), "\n")
  
  if ("parameters" %in% names(data_object)) {
    params <- data_object$parameters
    cat("Correlation structure:", params$correlation_structure, "\n")
    cat("Correlation parameter (rho):", params$rho, "\n")
    cat("Effect type:", params$effect_type, "\n")
    cat("Base effect size:", params$base_effect, "\n")
  }
  
  if ("true_ate" %in% names(data_object)) {
    cat("True average treatment effect:", round(data_object$true_ate, 4), "\n")
  }
  
  # Covariate correlation summary
  cor_matrix <- cor(data_object$X)
  diag(cor_matrix) <- NA
  cat("Covariate correlations:\n")
  cat("  Mean absolute correlation:", round(mean(abs(cor_matrix), na.rm = TRUE), 3), "\n")
  cat("  Max absolute correlation:", round(max(abs(cor_matrix), na.rm = TRUE), 3), "\n")
  cat("  Min absolute correlation:", round(min(abs(cor_matrix), na.rm = TRUE), 3), "\n")
  
  cat("==============================\n")
}