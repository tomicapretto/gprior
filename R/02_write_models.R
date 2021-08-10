library(rstan)
source(here::here("01_stan_code_atoms.R"))

write_model = function(program, path) {
  code = program$write()
  
  cat("Saving Stan code...\n")
  # Write .stan file
  con = file(paste0(path, ".stan"))
  writeLines(code, con)
  close(con)
  
  cat("Saving Stan compiled code...\n")
  # Compile and write model
  saveRDS(stan_model(model_code = code), paste0(path, ".rds"))
  cat("Done!\n")
}

# Model 1 -----------------------------------------------------------------
# g = n, log(sigma) ~ uniform.

program = StanProgram$new(
  Data$new(
    Number$new("n", "int"),
    Number$new("p", "int"),
    Number$new("g", "real"),
    Vector$new("y", "n"),
    Vector$new("mu_b", "p"),
    Matrix$new("mu_Sigma", "p", "p"),
    Matrix$new("X", "n", "p")
  ),
  Parameters$new(
    Vector$new("beta", "p"),
    Number$new("log_sigma_sq", "real")
  ),
  TransformedParameters$new(
    Assignment$new(
      Number$new("sigma", "real", lower=0),
      "sqrt(exp(log_sigma_sq))"
    )
  ),
  Model$new(
    SamplingStatement$new(
      "beta", "multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma)"
    ),
    SamplingStatement$new(
      "y", "normal(to_vector(X * beta), sigma)"
    )
  ),
  GeneratedQuantities$new(
    Vector$new("log_lik", "n"),
    ForLoop$new(
      "i", "1:n",
      Assignment$new("log_lik[i]", "normal_lpdf(y[i] | X[i] * beta, sigma)")
    )
  )
)

write_model(program, here::here("models", "model_1"))

# Model 2 -----------------------------------------------------------------
# g = n, sigma ~ exponential(1/sd(y))

program = StanProgram$new(
  Data$new(
    Number$new("n", "int"),
    Number$new("p", "int"),
    Number$new("g", "real"),
    Vector$new("y", "n"),
    Vector$new("mu_b", "p"),
    Matrix$new("mu_Sigma", "p", "p"),
    Matrix$new("X", "n", "p")
  ),
  Parameters$new(
    Vector$new("beta", "p"),
    Number$new("sigma", "real", lower = 0)
  ),
  Model$new(
    SamplingStatement$new(
      "beta", "multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma)"
    ),
    SamplingStatement$new(
      "y", "normal(to_vector(X * beta), sigma)"
    ),
    SamplingStatement$new(
      "sigma", "exponential(1 / sd(y))"
    )
  ),
  GeneratedQuantities$new(
    Vector$new("log_lik", "n"),
    Assignment$new(
      Number$new("sigma_p", "real"), 
      "exponential_rng(1 / sd(y))"
    ),
    Assignment$new(
      Vector$new("beta_p", "p"), 
      "multi_normal_rng(mu_b, g * pow(sigma_p, 2) * mu_Sigma)"
    ),
    Assignment$new(
      Array$new("y_p", "real", shape = "n"), 
      "normal_rng(to_vector(X * beta_p), sigma_p)"
    ),
    ForLoop$new(
      "i", "1:n",
      Assignment$new("log_lik[i]", "normal_lpdf(y[i] | X[i] * beta, sigma)")
    )
  )
)

write_model(program, here::here("models", "model_2"))

# Model 3 -----------------------------------------------------------------
# g = n, sigma ~ student_t(nu=4, sigma=sd(y))

program = StanProgram$new(
  Data$new(
    Number$new("n", "int"),
    Number$new("p", "int"),
    Number$new("g", "real"),
    Vector$new("y", "n"),
    Vector$new("mu_b", "p"),
    Matrix$new("mu_Sigma", "p", "p"),
    Matrix$new("X", "n", "p")
  ),
  Parameters$new(
    Vector$new("beta", "p"),
    Number$new("sigma", "real", lower = 0)
  ),
  Model$new(
    SamplingStatement$new(
      "beta", "multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma)"
    ),
    SamplingStatement$new(
      "y", "normal(to_vector(X * beta), sigma)"
    ),
    SamplingStatement$new(
      "sigma", "student_t(4, 0, sd(y))"
    )
  ),
  GeneratedQuantities$new(
    Vector$new("log_lik", "n"),
    Assignment$new(
      Number$new("sigma_p", "real"), 
      "fabs(student_t_rng(4, 0, sd(y)))"
    ),
    Assignment$new(
      Vector$new("beta_p", "p"), 
      "multi_normal_rng(mu_b, g * pow(sigma_p, 2) * mu_Sigma)"
    ),
    Assignment$new(
      Array$new("y_p", "real", shape = "n"), 
      "normal_rng(to_vector(X * beta_p), sigma_p)"
    ),
    ForLoop$new(
      "i", "1:n",
      Assignment$new("log_lik[i]", "normal_lpdf(y[i] | X[i] * beta, sigma)")
    )
  )
)

write_model(program, here::here("models", "model_3"))