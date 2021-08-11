library(here)
library(rstan)
source(here("02_prepare_simulation", "simulator.R"))
source(here("02_prepare_simulation", "data_generator.R"))

model = readRDS(here("models/model_1.rds"))

generator = DataGenerator$new(
  list("beta" = c(2, 0.8, -1.5, -0.3), "sigma" = 2),
  list("mu_Sigma" = TRUE, "g" = TRUE, "b_sd" = NULL, "rho" = 0.3)
)

SIZES = c(20, 40, 60)
REPS = 10

simulator = Simulator$new(model, generator)
simulator$make_plan(SIZES, REPS)
results = simulator$simulate()
