#' Generate Multivariate Normal Covariates
#'
#' Generates a matrix of correlated covariates from a multivariate normal distribution.
#'
#' @param n Sample size (number of observations)
#' @param p Number of covariates
#' @param rho Correlation coefficient between covariates (default: 0.3)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return An n x p matrix of covariates
#'
#' @examples
#' # Generate 100 observations with 10 covariates and correlation 0.5
#' X <- generate_covariates(n = 100, p = 10, rho = 0.5, seed = 123)
#' 
#' @export
#' @importFrom MASS mvrnorm
generate_covariates <- function(n = 100, p = 10, rho = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Input validation
  if (n <= 0 || p <= 0) stop("n and p must be positive integers")
  if (rho < -1 || rho > 1) stop("rho must be between -1 and 1")
  
  # Create correlation matrix
  mu <- rep(0, p)
  sigma <- matrix(rho, p, p)
  diag(sigma) <- 1
  
  # Generate data
  X <- MASS::mvrnorm(n, mu, sigma)
  
  # Standardize
  X <- scale(X)
  
  return(X)
}

#' Generate Response Variable
#'
#' Generates response variable with treatment effect and covariate relationships.
#'
#' @param X Covariate matrix (n x p)
#' @param treatment Treatment assignment vector (0/1)
#' @param beta Coefficient vector for covariates (default: rep(0.5, ncol(X)))
#' @param tau Treatment effect (default: 1)
#' @param sigma Error standard deviation (default: 1)
#' @param response_type Type of response function: "linear", "nonlinear1", "nonlinear2", "nonlinear3", "nonlinear4" (default: "linear")
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return Vector of response values
#'
#' @examples
#' X <- generate_covariates(100, 5, seed = 123)
#' treatment <- sample(c(0,1), 100, replace = TRUE)
#' y <- generate_response(X, treatment, tau = 1.5, seed = 123)
#'
#' @export
generate_response <- function(X, treatment, beta = NULL, tau = 1, sigma = 1, 
                              response_type = "linear", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Input validation
  if (length(treatment) != n) stop("Length of treatment must equal number of rows in X")
  if (!all(treatment %in% c(0, 1))) stop("Treatment must be binary (0/1)")
  
  # Default beta
  if (is.null(beta)) beta <- rep(0.5, p)
  if (length(beta) != p) stop("Length of beta must equal number of columns in X")
  
  # Calculate response based on type
  if (response_type == "linear") {
    g_fun <- X %*% beta
  } else if (response_type == "nonlinear1") {
    g_fun <- exp(X) %*% beta
  } else if (response_type == "nonlinear2") {
    g_fun <- exp(2 * X %*% beta / p)
  } else if (response_type == "nonlinear3") {
    g_fun <- sin(pi * X %*% beta / p)
  } else if (response_type == "nonlinear4") {
    g_fun <- (X %*% beta / p)^2
  } else {
    stop("Invalid response_type. Must be one of: linear, nonlinear1, nonlinear2, nonlinear3, nonlinear4")
  }
  
  # Generate response
  y <- as.vector(g_fun) + tau * treatment + sigma * rnorm(n)
  
  return(y)
}