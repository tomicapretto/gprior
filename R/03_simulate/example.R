library(here)
library(rstan)
source(here("02_prepare_simulation", "simulator.R"))
source(here("02_prepare_simulation", "data_generator.R"))
source(here("03_simulate", "settings.R"))

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
results = simulator$simulate()

sum(sapply(results, function(x) sum(x$sampler$time))) / 60 # ~ 18 minutes

