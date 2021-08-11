#' Obtain random samples from a Multivriate Normal Distribution with 
#' homogeneous correlation
#'
#' @param mu mean vector
#' @param sigma standard deviation vector
#' @param rho vector of length one with the common correlation
#' @param n number of samples to obtain
mvrnorm = function(mu, sigma, rho, n = 1) {
  # Omega: Correlation matrix
  # Sigma: Covariance matrix
  
  Omega = matrix(rho, length(mu), length(mu)) # p * p
  diag(Omega) = 1
  
  Sigma = diag(sigma) %*% Omega %*% diag(sigma) # p * p
  Chol = t(chol(Sigma)) # p * p, lower triangular (upper-triangular by default)
  
  univariate_samples = matrix(rnorm(n * length(mu)), ncol = n) # p * n
  samples = Chol %*% univariate_samples + mu # p * n
  
  return(t(samples)) # n * p
}

#' Generate data for particular simulation scenario
#'
#' @description
#' This class contains a function to generate samples and the parameters 
#' passed to that function. 
DataGenerator = R6::R6Class(
  "DataGenerator",
  public = list(
    params = NULL,
    mu_Sigma = NULL,
    g = NULL,
    b_sd = NULL,
    
    #' @field params A list with all the values of the parameters that determine
    #'  how the data are generated. This does not contain the sample size.
    initialize = function(params, settings) {
      self$params = params
      self$mu_Sigma = settings$mu_Sigma
      self$g = settings$g
      self$b_sd = settings$b_sd
      
      private$init(params, settings)
    },
    
    #' @description
    #' Obtain a random sample
    #' @param size The size of the random sample.
    data = function(size) {
      X = self$generate_X(size)
      y = self$generate_y(X, size)
      n = nrow(X)
      p = ncol(X)
      mu_b = rep(0, p)
      
      data = list("X" = X, "y" = y, "n" = n, "p" = p, "mu_b" = mu_b)
      
      if (self$mu_Sigma) {
        data[["mu_Sigma"]] = self$generate_Sigma(X)
      }
      
      if (self$g) {
        data[["g"]] = nrow(X)
      }
      
      if (!is.null(self$b_sd)) {
        data[["b_sd"]] = b_sd
      }
      return(data)
    },
    
    generate_X = function(size) {
      mvrnorm(private$x_mu, private$x_sigma, private$x_rho, size)
    },
    
    generate_y = function(X, size) {
      drop(X %*% self$params$beta + rnorm(size, sd = self$params$sigma))
    },
    
    generate_Sigma = function(X) {
      solve(t(X) %*% X)
    }
  ),
  
  private = list(
    x_mu = NULL,
    x_sigma = NULL,
    x_rho = NULL,
    init = function(params, settings) {
      private$x_mu = rep(0, length(params$beta))
      private$x_sigma = rep(1, length(params$beta))
      private$x_rho = settings$rho
    }
  )
)