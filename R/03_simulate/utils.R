simulate = function(MODEL, SCENARIO) {
  model = readRDS(here("models", MODEL$file))
  mu_Sigma = MODEL$mu_Sigma
  g = MODEL$g
  b_sd = MODEL$b_sd
  
  sigma = SCENARIO$SIGMA
  sizes = SCENARIO$SIZES
  reps = SCENARIO$REPS
  
  RESULT = list()
  
  for (betas in SCENARIO$BETAS) {
    for (rho in SCENARIO$RHO) {
      generator = DataGenerator$new(
        list("beta" = betas, "sigma" = sigma),
        list("mu_Sigma" = mu_Sigma, "g" = g, "b_sd" = b_sd, "rho" = rho)
      )
      
      simulator = Simulator$new(model, generator)
      simulator$make_plan(sizes, reps)
      
      result = simulator$simulate()
      RESULT[paste("beta", length(betas), "rho", rho, sep = "_")] = result
    }
  }
  return(results)
}



# for betas in BETAS (low, high)
# for rho in RHO
# 