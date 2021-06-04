Simulator = R6::R6Class(
  "Simulator",
  public = list(
    stan_model = NULL,
    pars = NULL,
    generate_data = NULL,
    probs = c(0.025, 0.975),
    results_params = list(),
    results_sampler = list(),
    
    initialize = function(stan_model, pars, generate_data) {
      self$stan_model = stan_model
      self$pars = pars
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
  )
)

library(rstan)
model = readRDS(here::here("models/model_1.rds"))

pars = c("beta", "sigma")
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


simulator = Simulator$new(model, pars, generate_data)
simulator$simulate(100)

# generate_data is a function that does not have any arguments
# because it already has all the information needed to generate the sample

# The 95% CI bounds are used to compute the length of the 95% intervals
# and wheter the true value is contained or not. Coverage should be close to 
# 95%.

# model@mode -> must be 0 to indicate it sampled.
