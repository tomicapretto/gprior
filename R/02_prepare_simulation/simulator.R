library(posterior)
library(rstan)
library(here)
source(here("utils", "utils.R"))

#' Simulator R6 class
#'
#' @description
#' This class wraps all the machinery to carry out simulations with a Stan
#' model and summarize results.
#'
#' @details
#' I'm playing with R6 documentation :')
Simulator = R6::R6Class(
  "Simulator",
  public = list(
    #' @field model Stan model.
    #' @field generator A DataGenerator object.
    #' @field params A list with names and values of the parameters in the model.
    #' @field sizes Sample sizes to use in the simulation.
    #' @field reps The number of repetitions.
    #' @field probs A vector of length two with the limits for the PI.
    model = NULL,
    generator = NULL,
    params = NULL,
    sizes = NULL,
    reps = NULL,
    messages = vector("character"),
    probs = c(0.025, 0.975),
   
    initialize = function(model, generator) {
      self$model = model
      self$generator = generator
      self$params = generator$params
    },
    
    #' @description
    #' Make a simulation plan. 
    #' @param sizes Sample sizes for the simulation.
    #' @param reps Repetitions for each sample size. Must be of the same size 
    #'   than \code{sizes}.      
    make_plan = function(sizes, reps) {
      self$sizes = sizes
      self$reps = reps
      private$planned = TRUE
    },
    
    #' @description 
    #' Simulate a single run.
    #' @param size Sample size used in this run.
    simulate_single = function(size) {
      # Generate data
      data = self$generator$data(size)
      
      # Run sampler
      fit = tryCatch({
        sampling(object = self$model, data = data, refresh = 0)
      }, 
      error = function(cnd) {
        cat("An error occurred while sampling!!!\n")
        cat(cnd$message, "\n")
        cnd$message
      })
      
      if (is.character(fit)) return(fit)
     
      # Return summary of the parameters and sampler diagnostics
      list(
        params = private$summarise_params(fit),
        sampler = get_sampler_params(fit, inc_warmup = FALSE),
        time = sum(get_elapsed_time(fit))
      )
    },
    
    #' @description
    #' Simulate \code{reps} repetitions for a sample of size \code{size}.
    #' @param size The sample size used in each of the repetitions.
    #' @param reps The number of repetitions.
    simulate_for_size = function(size) {
      first = TRUE
      print(glue::glue("Simulating with size: {size} and reps: {self$reps}\n"))
      pb = progress::progress_bar$new(total = self$reps)

      for (i in seq(self$reps)) {
        outcome = self$simulate_single(size)
        # When sampler fails, skip rest of loop and record message.
        if (is.character(outcome)) {
          self$messages = outcome
          next
        }
        if (first) {
          results = private$make_containers(outcome)
          first = FALSE
        }
        # Append parameter summaries
        for (param in names(results$params)) {
          results$params[[param]][i, ] = as.matrix(subset(outcome$params, variable == param)[-1])
        }
        
        # Append divergence count
        results$sampler$divergences[i] = private$get_divergences(outcome$sampler)
        results$sampler$time[i] = outcome$time
        pb$tick()
      }
      
      return(results)
    },
    
    #' @description
    #' Run simulation
    simulate = function() {
      if (!private$planned) stop("Plan the simulation first!!")
      results = lapply(self$sizes, self$simulate_for_size)
      names(results) = self$sizes
      results
    }
  ),
  private = list(
    planned = FALSE,
    
    #' @description
    #' Extract posterior samples from a fitted model.
    #' 
    #' This method only extracts the posterior samples for parameters in 
    #' \code{self$params} and combines them into a matrix.
    #' @param fit An object of class \code{stanfit}
    extract = function(fit) {
      x = do.call(cbind, extract(fit, pars = names(self$params)))
      colnames(x) = private$make_param_names()
      x
    },
    
    #' @description
    #' Make parameter names 
    #' 
    #' This method considers returns a character vector with names for the 
    #' parameters in the matrix with the posterior samples. It can also handle
    #' parameters with length larger than 1.
    make_param_names = function() {
      fun = function(name, value) {
        n = length(value)
        if (n > 1) {
          as.character(glue::glue("{name}[{seq_len(n)}]"))
        } else {
          name
        }
      }
      unname(unlist(mapply(fun, names(self$params), self$params)))
    },
    
    #' @description
    #' Get summaries for the parameters in the model.
    #' 
    #' @param fit An object of class \code{stanfit}.
    summarise_params = function(fit) {
      # Compute summaries using functions from the 'posterior' package.
      summary = private$summarise_draws(subset(as_draws(fit), names(self$params)))
      
      # Sort according to parameter names
      summary = summary[match(private$make_param_names(), summary$variable), ]
      
      # The following require the whole posterior and the true param values
      param_values = unname(unlist(self$params))
      posterior = private$extract(fit)
      
      # Compute RMSE
      rmse = private$get_rmse(posterior, param_values)
      # Compute P(beta_posterior < beta_true)
      prob = private$get_prob(posterior, param_values)
      # Compute bias 
      bias = private$get_bias(posterior, param_values)
      
      # Append computations
      summary$bias = bias
      summary$rmse = rmse
      summary$prob = prob
      
      # Derive coverage indicators
      summary$P50 = ifelse(between(param_values, summary$q25, summary$q75), 1, 0)
      summary$P80 = ifelse(between(param_values, summary$q10, summary$q90), 1, 0)
      summary$P90 = ifelse(between(param_values, summary$q5, summary$q95), 1, 0)
      summary$P95 = ifelse(between(param_values, summary$q2.5, summary$q97.5), 1, 0)
      
      drop = paste0("q", c(2.5, 5, 10, 25, 75, 90, 95, 97.5))
      summary = summary[!colnames(summary) %in% drop]
      summary
    },
    
    # TODO: What is this?
    get_coverage = function(draws, coef_true, probs) {
      # Probs is a list with vectors of length 2
      lapply(probs, function(x) {
        quantile()
      })
    },
    
    #' @description
    #' Compute RMSE (across posterior draws)
    #' 
    #' @param posterior A matrix with posterior samples. Each column represents
    #'   a coefficient.
    #' @param coef_true A vector with the true value of the coefficients in 
    #'   the model.
    get_rmse = function(posterior, coef_true) {
      sqrt(colMeans(sweep(posterior, 2, coef_true, function(x, y) (x - y) ^ 2)))
    },
    
    #' @description
    #' Compute P(beta_posterior < beta_true)
    #' 
    #' @param posterior A matrix with posterior samples. Each column represents
    #'   a coefficient.
    #' @param coef_true A vector with the true value of the coefficients in 
    #'   the model.
    get_prob = function(posterior, coef_true) {
      # I dont remember the name of this quantity
      colMeans(sweep(posterior, 2, coef_true, FUN = "<"))
    },
    
    #' @description
    #' Compute Bias (of the posterior mean)
    #' 
    #' @param posterior A matrix with posterior samples. Each column represents
    #'   a coefficient.
    #' @param coef_true A vector with the true value of the coefficients in 
    #'   the model.
    get_bias = function(posterior, coef_true) {
      colMeans(sweep(posterior, 2, coef_true, FUN = "-"))
    },
    
    
    #' @description
    #' Get summaries from posterior draws.
    #' 
    #' This method wraps \code{summarise_draws} from \pkg{posterior}.
    #' @param draws An object of class \code{draws} from \pkg{posterior}.
    summarise_draws = function(draws) {
      probs = c(0.025, 0.05, 0.1,  0.25, 0.75, 0.9, 0.95, 0.975)
      summarise_draws(draws, mean, sd, ess_bulk, rhat, ~quantile2(.x, probs))
    },
    
    #' @description
    #' Count the number of divergences in each simulation run.
    #'
    #' @param sampler_results The outcome of \code{get_sampler_params(fit)}.
    get_divergences = function(sampler_results) {
      sum(sapply(sampler_results, function(x) sum(x[, "divergent__"])))
    },
    
    #' @description
    #' Create objects that will contain the simulation results.
    #'
    #' @param results The outcome of a single simulation.
    make_containers = function(results) {
      # How many parameters are in the model
      n_params = nrow(results$params)
      
      # How many columns we put in the result matrix for each parameter.
      n_col_params = ncol(results$params) - 1
      
      # How many rows we put in the results matrix for each parameter.
      n_row_params = self$reps
      
      # Container for the summaries for the marginal posteriors
      params = replicate(
        n_params, 
        matrix(
          nrow = n_row_params, 
          ncol = n_col_params,
          dimnames = list(NULL, colnames(results$params)[-1])
        ), 
        simplify = FALSE
      )
      names(params) = results$params$variable
      
      # Information about the sampling process.
      sampler = list(
        divergences = vector("numeric", self$reps),
        time = vector("numeric", self$reps)
      )
      
      list("params" = params, "sampler" = sampler)
    }
  )
)

# model@mode -> must be 0 to indicate it sampled.

## Things to add
# yhat -> mu = X %*% beta, where beta is a (vector) draw from the posterior
#         sd = The one obtained from the posterior.
#         rnorm(1, mu, sd)

# Cross validation
# See https://discourse.mc-stan.org/t/bayesian-rmse/13293/4