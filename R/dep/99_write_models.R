library(rstan)
source(here::here("01_stan_model_code.R"))

write_file = function(content, path) {
  con = file(path)
  writeLines(content, con)
  close(con)
}

# Model 1 -----------------------------------------------------------------
# g = n, log(sigma) ~ uniform.

# Create model object
stan_code = StanModelCode$new()

# Data
stan_code$add_lines(
  block = "data", 
  lines = list(
    list(name = "n", type = "int"), 
    list(name = "p", type = "int"), 
    list(name = "g", type = "real"),
    list(name = "y", type = "vector", rows = "n"),
    list(name = "mu_b", type = "vector", rows = "p"),
    list(name = "mu_Sigma", type = "matrix", rows = "p", cols = "p"),
    list(name = "X", type = "matrix", "rows" = "n", "cols" = "p")
  )
)

# Parameters
stan_code$add_lines(
  block = "parameters",
  lines = list(
    list(name = "beta", type = "vector", rows = "p"),
    list(name = "log_sigma_sq", type = "real")
  )
)

# Transformed parameters
stan_code$add_line(
  block = "transformed parameters",
  line = list(
    name = "sigma", 
    type = "real", 
    lower = 0, 
    assignment = "sqrt(exp(log_sigma_sq))"
  )
)

# Model
stan_code$add_lines(
  block = "model",
  lines = list(
    list(
      name = "beta", 
      type = "dist", 
      dist = "multi_normal", 
      args = c("mu_b", "g * pow(sigma, 2) * mu_Sigma")
    ),
    list(
      name = "y", 
      type = "dist", 
      dist = "normal", 
      args = c("to_vector(X * beta)", "sigma")
    )
  )
)

# Log likelihood, for loo computation
stan_code$add_lines(
  block = "generated quantities",
  lines = list(
    list(name = "log_lik", type = "vector", rows = "n"),
    list(
      name = "log_lik_loop",
      type = "loop", 
      var = "i",
      seq = "1:n",
      expr = list(
        name = "log_lik[i]",
        props = list(
          type = "assignment",
          expr = "normal_lpdf(y[i] | X[i] * beta, sigma)"
        )
      )
    )
  )
)

# Obtain code, write stan file and the compiled model
model_1_code = stan_code$make_stan_code()
write_file(model_1_code, here::here("models", "model_1.stan"))
saveRDS(stan_model(model_code = model_1_code), here::here("models", "model_1.rds"))

# Model 2 -----------------------------------------------------------------
# g = n, sigma ~ exponential(1/sd(y))

stan_code$remove_block("transformed parameters")
stan_code$remove_line("parameters", "log_sigma_sq")
stan_code$remove_line("model", "log_sigma_sq")
stan_code$remove_line("generated quantities", "log_sigma_sq_p")

stan_code$add_line(
  "parameters", 
  list(
    name = "sigma", 
    type = "real", 
    lower = 0
  )
)
stan_code$add_line(
  "model", 
  list(
    name = "sigma", 
    type = "dist", 
    dist = "exponential", 
    args = "1 / sd(y)"
  )
)

# Generated quantities
stan_code$add_lines(
  block = "generated quantities",
  lines = list(
    list(
      name = "sigma_p", 
      type = "real", 
      assignment = "exponential_rng(1 / sd(y))"
    ),
    list(
      name = "beta_p", 
      type = "vector", 
      rows = "p", 
      assignment = "multi_normal_rng(mu_b, g * pow(sigma_p, 2) * mu_Sigma)"
    ),
    list(
      name = "y_p", 
      type = "real", 
      dims = "[n]",
      assignment = "normal_rng(to_vector(X * beta_p), sigma_p)"
    )
  )
)

model_2_code = stan_code$make_stan_code()
write_file(model_2_code, here::here("models", "model_2.stan"))
saveRDS(stan_model(model_code = model_2_code), here::here("models", "model_2.rds"))

# Model 3 -----------------------------------------------------------------
# g = n, sigma ~ student_t(nu=4, sigma=sd(y))
stan_code$add_line(
  "model", 
  list(
    name = "sigma", 
    type = "dist",
    dist = "student_t", 
    args = c("4", "0", "sd(y)")
  )
)
stan_code$add_line(
  "generated quantities", 
  list(
    name = "sigma_p", 
    type = "real", 
    assignment = "fabs(student_t_rng(4, 0, sd(y)))"
  )
)
model_3_code = stan_code$make_stan_code()
write_file(model_3_code, here::here("models", "model_3.stan"))
saveRDS(stan_model(model_code = model_3_code), here::here("models", "model_3.rds"))

# Model 4 -----------------------------------------------------------------
# g ~ student_t(3, 0, 4), log(sigma) ~ uniform
stan_code$remove_line("data", "g")
stan_code$remove_line("model", "sigma")
stan_code$remove_line("parameters", "sigma")
stan_code$remove_line("generated quantities", "sigma_p")

stan_code$add_lines(
  block = "parameters",
  lines = list(
    list(name = "log_sigma_sq", type = "real"),
    list(name = "g", type = "real", lower = 0)
  )
)

# Add transformed parameters
stan_code$add_line(
  block = "transformed parameters",
  line = list(name = "sigma", type = "real", lower = 0, assignment = "sqrt(exp(log_sigma_sq))")
)

# Add prior for g
stan_code$add_line(
  "model", 
  list(name = "g", type = "dist", dist = "student_t", args = c(3, 0, 4))
)

# Generated quantities
stan_code$remove_block("generated quantities")
stan_code$sort_blocks()

model_4_code = stan_code$make_stan_code()
write_file(model_4_code, here::here("models", "model_4.stan"))
saveRDS(stan_model(model_code = model_4_code), here::here("models", "model_4.rds"))

# Model 5 -----------------------------------------------------------------
# g ~ student_t(3, 0, 4), sigma ~ exp(1 / sd(y))
stan_code$remove_block("transformed parameters")
stan_code$remove_line("parameters", "log_sigma_sq")
stan_code$remove_line("model", "log_sigma_sq")

stan_code$add_line("parameters", list(name = "sigma", type = "real", lower = 0))
stan_code$add_line("model", list(name = "sigma", type = "dist", dist = "exponential", args = "1 / sd(y)"))

stan_code$remove_line("generated quantities", "log_sigma_sq_p")

stan_code$add_lines(
  block = "generated quantities",
  lines = list(
    list(
      name = "sigma_p", 
      type = "real", 
      lower = 0,
      assignment = "exponential_rng(1/sd(y))"
    ),
    list(
      name = "g_p", 
      type = "real", 
      assignment = "fabs(student_t_rng(3, 0, 4))"
    ),
    list(
      name = "beta_p", 
      type = "vector", 
      rows = "p", 
      assignment = "multi_normal_rng(mu_b, g_p * pow(sigma_p, 2) * mu_Sigma)"
    ),
    list(
      name = "y_p", 
      type = "real", 
      dims = "[n]",
      assignment = "normal_rng(to_vector(X * beta_p), sigma_p)"
    )
  )
)

model_5_code = stan_code$make_stan_code()
write_file(model_5_code, here::here("models", "model_5.stan"))
saveRDS(stan_model(model_code = model_5_code), here::here("models", "model_5.rds"))

# Model 6 -----------------------------------------------------------4------
# g ~ student_t(3, 0, 4), sigma ~ student_t(nu=4, sigma=sd(y))
stan_code$add_line(
  "model", 
  list(name = "sigma", type = "dist", dist = "student_t", args = c("4", "0", "sd(y)"))
)
stan_code$add_line(
  "generated quantities", 
  list(
    name = "sigma_p", 
    type = "real", 
    assignment = "fabs(student_t_rng(4, 0, sd(y)))"
  )
)

model_6_code = stan_code$make_stan_code()
write_file(model_6_code, here::here("models", "model_6.stan"))
saveRDS(stan_model(model_code = model_6_code), here::here("models", "model_6.rds"))

# Model 7 -----------------------------------------------------------------

# This aims to use independent normal priors with sd = b_sd * sd(y), where 'b_sd'
# is a constant we specify. b_sd = 1 is suggested by Gelman and c = 2.5 is used
# in rstanarm. Note that rstanarm docs use sd(y)/sd(x). We don't need to divide
# by "x" because the predictors are standardized.
stan_code$remove_line("data", "mu_Sigma")
stan_code$remove_line("parameters", "g")
stan_code$remove_line("model", "g")
stan_code$add_line(
  "data",
  list(
    name = "b_sd",
    type = "real",
    lower = 0
  )
)

stan_code$add_lines(
  "model", 
  list(
    list(
      name = "beta", 
      type = "dist", 
      dist = "normal", 
      args = c("mu_b", "b_sd * sd(y)")
    ),
    list(
      name = "sigma",
      type = "dist", 
      dist = "exponential", 
      args = "1 / sd(y)"
    )
  )
)

stan_code$remove_line("generated quantities", "g_p")
stan_code$add_lines(
  "generated quantities", 
  list(
    list(
      name = "sigma_p", 
      type = "real",
      lower = 0,
      assignment = "exponential_rng(1 / sd(y))"
    ),
    list(
      name = "beta_p", 
      type = "real", 
      dims = "[p]",
      assignment = "normal_rng(mu_b, b_sd * sd(y))"
    ),
    list(
      name = "y_p", 
      type = "real", 
      dims = "[n]",
      assignment = "normal_rng(to_vector(X * to_vector(beta_p)), sigma_p)"
    )
  )
)

model_7_code = stan_code$make_stan_code()
write_file(model_7_code, here::here("models", "model_7.stan"))
saveRDS(stan_model(model_code = model_7_code), here::here("models", "model_7.rds"))

# Model 8 -----------------------------------------------------------------
# Flat priors everywhere!

stan_code$remove_block("generated quantities")
stan_code$remove_line("data", "b_sd")
stan_code$remove_line("model", "beta")
stan_code$remove_line("model", "sigma")

model_8_code = stan_code$make_stan_code()
write_file(model_8_code, here::here("models", "model_8.stan"))
saveRDS(stan_model(model_code = model_8_code), here::here("models", "model_8.rds"))


