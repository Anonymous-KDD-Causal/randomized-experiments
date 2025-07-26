#' Data Generation Functions for Experimental Design
#' 
#' Functions to generate synthetic data for testing experimental design methods.
#' 
#' @author Your Name
#' @version 1.0.0

#' Generate Correlated Covariates
#'
#' @param d Number of covariates
#' @param rho Correlation coefficient
#' @param n Sample size (default: 100)
#' @param seed Random seed (default: NULL)
#' @return Matrix of covariates (n x d)
#' @export
#' @examples
#' # Generate 100 observations with 5 correlated covariates
#' X <- generate_covariates(d = 5, rho = 0.3, n = 100, seed = 123)
generate_covariates <- function(d, rho, n = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create correlation matrix
  mu <- rep(0, d)
  sigma <- matrix(rho, d, d)
  diag(sigma) <- 1
  
  # Generate multivariate normal data
  X <- MASS::mvrnorm(n, mu, sigma)
  return(X)
}

#' Generate Response Variable
#'
#' @param X Covariate matrix (n x p)
#' @param w Treatment allocation vector (0/1)
#' @param beta Coefficient vector for covariates
#' @param tau Treatment effect (default: 1)
#' @param std_err Error standard deviation (default: 1)
#' @param type Response type: 'linear', 'nonlinear1', 'nonlinear2', 'nonlinear3', 'nonlinear4'
#' @param seed Random seed (default: NULL)
#' @return Response vector
#' @export
#' @examples
#' # Generate linear response
#' X <- generate_covariates(5, 0.3, 100)
#' w <- sample(c(0, 1), 100, replace = TRUE)
#' y <- generate_response(X, w, beta = rep(0.5, 5), tau = 1, type = 'linear')
generate_response <- function(X, w, beta, tau = 1, std_err = 1, type = 'linear', seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Calculate systematic component based on type
  if (type == 'linear') {
    g_fun <- X %*% beta
  } else if (type == 'nonlinear1') {
    g_fun <- exp(X) %*% beta
  } else if (type == 'nonlinear2') {
    g_fun <- exp(2 * X %*% beta / d)
  } else if (type == 'nonlinear3') {
    g_fun <- sin(pi * X %*% beta / d)
  } else if (type == 'nonlinear4') {
    g_fun <- (X %*% beta / d)^2
  } else {
    stop("Invalid response type. Choose from: 'linear', 'nonlinear1', 'nonlinear2', 'nonlinear3', 'nonlinear4'")
  }
  
  # Generate response
  y <- g_fun + tau * w + std_err * rnorm(n)
  return(as.numeric(y))
}

#' Generate Complete Experimental Dataset
#'
#' @param n Sample size
#' @param p Number of covariates
#' @param rho Correlation between covariates
#' @param tau Treatment effect
#' @param response_type Type of response function
#' @param seed Random seed
#' @return List containing X (covariates), w (allocation), y (response)
#' @export
#' @examples
#' # Generate complete dataset
#' data <- generate_experiment_data(n = 100, p = 5, rho = 0.3, tau = 1)
generate_experiment_data <- function(n = 100, p = 5, rho = 0.3, tau = 1, 
                                     response_type = 'linear', seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  X <- generate_covariates(d = p, rho = rho, n = n)
  X <- scale(X)  # Standardize
  
  # Generate balanced allocation
  w <- sample(c(rep(0, n/2), rep(1, n/2)))
  
  # Generate response
  beta <- rep(0.5, p)  # Default coefficient vector
  y <- generate_response(X, w, beta, tau, type = response_type)
  
  return(list(
    X = X,
    w = w,
    y = y,
    n = n,
    p = p,
    rho = rho,
    tau = tau,
    response_type = response_type
  ))
}
