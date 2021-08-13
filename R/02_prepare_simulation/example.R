library(here)
library(rstan)
source(here("02_prepare_simulation", "simulator.R"))
source(here("02_prepare_simulation", "data_generator.R"))

# Generate values for the coefficients ------------------------------------
set.seed(121195)

get_betas_for_size = function(n, SDs = c(2, 4)) {
  list("low" = rnorm(n, sd = SDs[1]), "high" = rnorm(n, sd = SDs[2]))
}

BETAS_SIZES =  c(3, 8, 16, 20, 30, 36, 50, 70, 100, 150, 180, 400)
BETAS = lapply(BETAS_SIZES, get_betas_for_size)
names(BETAS) = BETAS_SIZES


# Define other misc parameters --------------------------------------------

# Homogeneous correlation
RHO = c(0.1, 0.3, 0.6, 0.9)

# Sigma
SIGMA = 2

SIZES = c(20, 40, 80, 120, 500, 2000)

# Number of repetitions
REPS = 200

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

MODEL = MODELS[[1]]
model = readRDS(here("models", MODEL$file))

SCENARIO = SCENARIOS[[1]]
BETA = SCENARIO$BETAS[[1]]$low
SIGMA = SCENARIO$SIGMA
RHO = SCENARIO$RHO[[1]]
SIZES = SCENARIO$SIZES
REPS = SCENARIO$REPS

generator = DataGenerator$new(
  list("beta" = BETA, "sigma" = SIGMA),
  list(
    "mu_Sigma" = MODEL$mu_Sigma, 
    "g" = MODEL$g, 
    "b_sd" = MODEL$b_sd, 
    "rho" = RHO
  )
)

simulator = Simulator$new(model, generator)
simulator$make_plan(SIZES, REPS)
results_1 = simulator$simulate()

results_1$`20`$sampler$time
  
sum(sapply(results_1, function(x) sum(x$sampler$time))) / 60
# ~ 18 minutes
