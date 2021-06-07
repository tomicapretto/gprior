Simulator = R6::R6Class(
  "Simulator",
  public = list(
    stan_model = NULL,
    pars = NULL,
    parsv = NULL,
    generate_data = NULL,
    probs = c(0.025, 0.975),
    results_params = list(),
    results_sampler = list(),
    
    initialize = function(stan_model, pars, parsv, generate_data) {
      # stan_model: A stan compiled model
      # pars: Names of the parameters in the model (e.g. beta and sigma)
      # parsv: Named list with the true values of the parameters in the model
      # generate_data: Function without arguments that generates data passed to
      #                rstan::sampling()
      self$stan_model = stan_model
      self$pars = pars
      self$parsv = parsv
      self$generate_data = generate_data
    },
    
    simulate = function(reps) {
      first = TRUE
      pb = progress::progress_bar$new(total = reps)
      for (i in seq(reps)) {
        result = self$simulate_single()
        if (first) {
          self$make_containers(result, reps)
          first = FALSE
        }
        # Append parameter summaries
        for (param in names(self$results_params)) {
          self$results_params[[param]][i, ] = result$params[param, ]
        }
        # Append divergence count
        self$results_sampler$divergences[i] = sum(sapply(
          result$sampler, 
          function(x) sum(x[, "divergent__"])
        ))
        pb$tick()
      }
      
      # Append other summaries such as bias, mse, coverage.
      for (par in names(self$parsv)) {
        value = self$parsv[[par]]
        if (length(value) > 1) {
          for (i in seq(length(value))) {
            name = glue::glue("{par}[{i}]")
            self$results_params[[name]] = private$enrich_summary(
              self$results_params[[name]], 
              value[[i]]
            )
          }
        } else {
          self$results_params[[par]] = private$enrich_summary(
            self$results_params[[par]], 
            value
          )
        }
      }
    },
  
    simulate_single = function() {
      # Generate data
      data = self$generate_data()
      # Run sampler
      fit = sampling(object = self$stan_model, data = data, refresh = 0)
      # Return summary of the parameters and sampler diagnostics
      return( 
        list(
          params = summary(fit, pars = self$pars, probs = self$probs)$summary,
          sampler = get_sampler_params(fit, inc_warmup = FALSE)
        )
      )
    },
    
    make_containers = function(results, reps) {
      n_params = nrow(results$params)
      n_col_params = ncol(results$params)
      n_row_params = reps
      
      # Summaries about the marginal posteriors
      self$results_params = replicate(
        n_params, 
        matrix(
          nrow = n_row_params, 
          ncol = n_col_params,
          dimnames = list(NULL, colnames(results$params))
        ), 
        simplify = FALSE
      )
      names(self$results_params) = rownames(results$params)
      
      # Information about the sampling process.
      # For now, only keep the count of divergences
      self$results_sampler = list(
        divergences = vector("numeric", reps)
      )
    }
  ),
  active = list(
    results = function() {
      list("params" = self$results_params, "sampler" = self$results_sampler)
    }
  ),
  private = list(
    # Take a summary matrix, adds bias, mse and coverage (0 or 1 for the PI)
    enrich_summary = function(mm, value) {
      bounds = paste0((self$probs * 100), "%")
      bias = mm[, "mean"] - value
      mse = bias ^ 2 + mm[, "sd"] ^ 2
      coverage = ifelse(value >= mm[, bounds[1]] | value <= mm[, bounds[2]], 1, 0)
      cbind(mm, bias, mse, coverage)
    }
  )
)

library(rstan)
model = readRDS(here::here("models/model_1.rds"))

generate_data = function() {
  SIZE = 30

  x1 = rnorm(SIZE)
  x2 = 0.35 + 0.15 * x1 + rnorm(SIZE)
  x3 = 0.5 + 0.3 * x1 + 0.2 * x2 + rnorm(SIZE)
  x4 = -0.7 + 0.2 * x1 + 0.1 * x3 + rnorm(SIZE)
  
  X = cbind(x1, x2, x3, x4)
  X = scale(X)
  b_true = c(2, 0.8, -1.5, -0.3)
  sigma_true = 2
  y = X %*% b_true + rnorm(SIZE, sd=sigma_true)
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

pars = c("beta", "sigma")
parsv = list("beta" = c(2, 0.8, -1.5, -0.3), "sigma" = 2)
simulator = Simulator$new(model, pars, parsv, generate_data)
simulator$simulate(50)

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