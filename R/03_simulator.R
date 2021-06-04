Simulator = R6::R6Class(
  "Simulator",
  public = list(
    stan_model = NULL,
    pars = NULL,
    data_f = NULL,
    data_args = NULL,
    probs = c(0.025, 0.975),
    
    initialize = function(stan_model, pars, data_f, data_args) {
      self$stan_model = stan_model
      self$pars = pars
      self$data_f = data_f
      self$data_args = data_args
    },
    
    simulate = function(reps) {
      first = TRUE
      for (i in seq(reps)) {
        result = self$simulate_single()
        if (first) {
          # Use the first run to allocate space for the parameter values
          self$make_containers(result, reps)
          first = FALSE
        }
        self$append_result(result)
      }
      
      
    },
    
    simulate_single = function() {
      # Generate data
      data = do.call(self$data_f, self$data_args)
      # Run sampler
      fit = sampling(object = self$stan_model, data = data, referesh = 0)
      # Compute summary
      summary(fit, pars = self$pars, probs = self$probs)$summary
      # Obtain sampler params 
    },
    
    make_containers = function(results, reps) {
      n_params = nrow(results$summary)
      n_col_params = ncol(results$summary)
      n_row_params = reps * 4
      
      # Summaries about the marginal posteriors
      self$results_params = replicate(
        n_params, 
        matrix(nrow = n_row_params, ncol = n_col_params), 
        simplify = FALSE
      )
      names(self$results_params) = rownames(results$summary)
      
      # Information about the sampling process.
      self$results_sampler = matrix(nrow = n_row_params * 4, ncol = 2)
      
      
    },
    
    get_results = function() {
      
    }
    
  ),
  
  private = list(
    
  )
)


# Things to consider when working with the output of stan models


# Each column is a parameter, rows are all the draws (stacked)
# matrix_of_draws <- as.matrix(fit)

# Obtain a summary of the fit
# fit_summary <- summary(fit)

# fit_summary$summary has info for chains merged.
# Here we are going to take
# mean, sd, 95% CI bounds, n_eff, and Rhat.

# The 95% CI bounds are used to compute the length of the 95% intervals
# and wheter the true value is contained or not. Coverage should be close to 
# 95%.

# We get a list with one matrix per chain (info is not combined)
# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# sampler_params_chain1 <- sampler_params[[1]]

# Interesting things to compute
# Sum of divergences using divergent__
# Pct of accepted proposals? Don't know if this is important.
# Ask Paul about other quantities.

# model@mode -> must be 0 to indicate it sampled.

summary(fit, pars = c("beta", "sigma"), probs = c(0.025, 0.975))$summary

sampler_params = get_sampler_params(fit, inc_warmup = FALSE)
