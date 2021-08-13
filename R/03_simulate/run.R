library(here)
library(rstan)
source(here("02_prepare_simulation", "simulator.R"))
source(here("02_prepare_simulation", "data_generator.R"))
source(here("03_simulate", "settings.R"))
source(here("03_simulate", "utils.R"))


MODEL = MODELS[[1]]
SCENARIO = SCENARIOS[[1]] 

start <- Sys.time()

RESULTS = list()
for (i in seq_along(MODELS)) {
  cat("Simulating for Model", i, "\n")
  RESULTS_MODEL = list()
  for (j in seq_along(SCENARIOS)) {
    cat("Simulating for Scenario", j, "\n")
    RESULTS_MODEL[[paste0("scenario_", j)]] = simulate(MODEL, SCENARIO)
  }
  RESULTS[[paste0("model_", i)]] = RESULTS_MODEL
}

end <- Sys.time()

end - start
