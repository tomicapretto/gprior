library(rstan)
library(posterior)

#' R6 Class Representing a Simulator
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
    #' Make simulating plan. 
    #' @param sizes Sample sizes to use in the simulation.
    #' @param reps Repetitions used with each sample size. Must be of the same size than sizes.      
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
    #' Simulate 'reps' repetitions for a sample of size 'size'.
    #' @param size The size used in each of the repetitions.
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
    
    # Main function
    simulate = function() {
      if (!private$planned) stop("Plan the simulation first!!")
      results = lapply(self$sizes, self$simulate_for_size)
      names(results) = self$sizes
      results
    }
  ),
  private = list(
    planned = FALSE,
    
    summarise_params = function(fit) {
      # First of all, compute summaries using 'posterior' functions
      draws = subset(as_draws(fit), names(self$params))
      summary = private$summarise_draws(draws)
      
      # Compute RMSE
      rmse = private$make_rmse_vector(fit, draws)
      
      # Sort RMSE according to order in summary and append to tibble
      summary$rmse = rmse[match(summary$variable, names(rmse))]
      summary
    },
    
    make_rmse_vector = function(fit, draws) {
      # Create empty vector
      rmse = vector("numeric", length(unlist(self$params)))
      # Get names of the empty vector
      new_names = sapply(names(self$params), function(name) {
        times = length(self$params[[name]])
        if (times > 1) {
          as.character(glue::glue("{name}[{seq_len(times)}]"))
        } else {
          name
        }
      })
      names(rmse) = unlist(unname(new_names))
      # Compute RMSE
      start = 1
      end = 0
      for (name in names(self$params)) {
        coef_true = self$params[[name]]
        end = end + length(coef_true)
        draws = extract(fit, name)[[name]]
        draws = matrix(draws, dim(draws)[1])
        rmse[start:end] = private$get_rmse(draws, coef_true)
        start = end + 1
      }
      rmse
    },
    
    get_rmse = function(draws, coef_true) {
      sqrt(colMeans(sweep(draws, 2, coef_true, function(x, y) (x - y) ^ 2)))
    },

    summarise_draws = function(draws) {
      probs = c(0.025, 0.05, 0.1, 0.2, 0.25, 0.75, 0.8, 0.9, 0.95, 0.975)
      summarise_draws(draws, mean, sd, ess_bulk, rhat, ~quantile2(.x, probs))
    },
    
    # Take a summary matrix, adds bias, mse and coverage (0 or 1 for the PI)
    enrich_summary = function(mm, value) {
      bounds = paste0((self$probs * 100), "%")
      bias = mm[, "mean"] - value
      mse = bias ^ 2 + mm[, "sd"] ^ 2
      coverage = ifelse(value >= mm[, bounds[1]] | value <= mm[, bounds[2]], 1, 0)
      cbind(mm, bias, mse, coverage)
    },
    
    #' @description
    #' Takes the outcome of `get_sampler_params(fit)` and returns the count of
    #' divergences in the sampling process.
    get_divergences = function(sampler_results) {
      sum(sapply(sampler_results, function(x) sum(x[, "divergent__"])))
    },
    
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

library(rstan)
model = readRDS(here::here("models/model_1.rds"))

generate_data = function(size, params) {
  
  # This has to be a call to mvnorm
  x1 = rnorm(size)
  x2 = 0.35 + 0.15 * x1 + rnorm(size)
  x3 = 0.5 + 0.3 * x1 + 0.2 * x2 + rnorm(size)
  x4 = -0.7 + 0.2 * x1 + 0.1 * x3 + rnorm(size)
  
  X = cbind(x1, x2, x3, x4)
  X = scale(X)
  y = X %*% params$beta + rnorm(size, sd=params$sigma)
  g = nrow(X)
  Sigma = solve(t(X) %*% X)
  
  list(
    n = nrow(X),
    p = ncol(X),
    g = nrow(X),
    y = as.vector(y),
    X = X,
    mu_Sigma = Sigma,
    mu_b = rep(0, 4)
  )
}

DataGenerator = R6::R6Class(
  "Simulator",
  public = list(
    make_data = NULL,
    params = NULL,
    initialize = function(make_data, params) {
      self$make_data = make_data
      self$params = params
    },
    data = function(size) {
      self$make_data(size, self$params)
    }
  )
)

params = list("beta" = c(2, 0.8, -1.5, -0.3), "sigma" = 2)
generator = DataGenerator$new(generate_data, params)

SIZES = c(20, 40, 60)
REPS = 10

simulator = Simulator$new(model, generator)
simulator$make_plan(SIZES, REPS)
results = simulator$simulate()



# See `simulator$results`

# model@mode -> must be 0 to indicate it sampled.

## Notes
# * Bias: is actually the bias of the posterior mean.
# * MSE: same than above
# * Coverage: It is 1 if the true value of the param is contained in the 95% PI.
# * yhat?
# * RMSE?

## Things to add
# yhat -> mu = X %*% beta, where beta is a (vector) draw from the posterior
#         sd = The one obtained from the posterior.
#         rnorm(1, mu, sd)
# RMSE -> sqrt(mean((yhat - y) ^ 2)) for each posterior sample
#         Then I have a posterior of the RMSE, same for MSE, etc

# Cross validation

# See https://discourse.mc-stan.org/t/bayesian-rmse/13293/4