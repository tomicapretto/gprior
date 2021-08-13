# Generate values for the regression coefficients -------------------------
set.seed(121195)

# get_betas_for_size = function(n, SDs = c(2, 4)) {
#   list("low" = rnorm(n, sd = SDs[1]), "high" = rnorm(n, sd = SDs[2]))
# }

BETAS_SIZES = c(3, 8, 16, 20, 30, 36, 50, 70, 100, 150, 180, 400)
BETAS = lapply(BETAS_SIZES, rnorm, sd = 2)
names(BETAS) = BETAS_SIZES

# Define other misc parameters --------------------------------------------
# Homogeneous correlation
RHO = c(0.1, 0.3, 0.6, 0.9)

# Sigma
SIGMA = 2

# Sample sizes
SIZES = c(20, 40, 80, 120, 500, 2000)

# Number of repetitions
REPS = 4

# Define settings ---------------------------------------------------------
SCENARIOS = list(
  list(
    "SIZES" = SIZES,
    "BETAS" = BETAS[1:3],
    "SIGMA" = SIGMA,
    "RHO" = RHO,
    "REPS" = REPS
  ),
  list(
    "SIZE" = SIZES[2:6],
    "BETAS" = BETAS[1:6],
    "SIGMA" = SIGMA,
    "RHO" = RHO,
    "REPS" = REPS
  ),
  list(
    "SIZE" = SIZES[2:6],
    "BETAS" = BETAS[1:8],
    "SIGMA" = SIGMA,
    "RHO" = RHO,
    "REPS" = REPS
  ),
  list(
    "SIZE" = SIZES[3:6],
    "BETAS" = BETAS[1:9],
    "SIGMA" = SIGMA,
    "RHO" = RHO,
    "REPS" = REPS
  ),
  list(
    "SIZE" = SIZES[4:6],
    "BETAS" = BETAS,
    "SIGMA" = SIGMA,
    "RHO" = RHO,
    "REPS" = REPS
  )
)

MODELS = list(
  list(
    "file" = "model_1.rds",
    "mu_Sigma" = TRUE,
    "g" = TRUE,
    "b_sd" = NULL
  ),
  list(
    "file" = "model_2.rds",
    "mu_Sigma" = TRUE,
    "g" = TRUE,
    "b_sd" = NULL
  ),
  list(
    "file" = "model_3.rds",
    "mu_Sigma" = TRUE,
    "g" = TRUE,
    "b_sd" = NULL
  ),
  list(
    "file" = "model_4.rds",
    "mu_Sigma" = TRUE,
    "g" = FALSE,
    "b_sd" = NULL
  ),
  list(
    "file" = "model_5.rds",
    "mu_Sigma" = TRUE,
    "g" = FALSE,
    "b_sd" = NULL
  ),
  list(
    "file" = "model_6.rds",
    "mu_Sigma" = TRUE,
    "g" = FALSE,
    "b_sd" = NULL
  ),
  list(
    "file" = "model_7.rds",
    "mu_Sigma" = FALSE,
    "g" = FALSE,
    "b_sd" = 1
  ),
  list(
    "file" = "model_7.rds",
    "mu_Sigma" = FALSE,
    "g" = FALSE,
    "b_sd" = 2.5
  ),
  list(
    "file" = "model_8.rds",
    "mu_Sigma" = FALSE,
    "g" = FALSE,
    "b_sd" = NULL
  )
)
