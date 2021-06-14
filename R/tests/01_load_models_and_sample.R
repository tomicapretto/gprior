library(rstan)
library(here)

# Toy sample --------------------------------------------------------------
SIZE = 30
set.seed(1995)

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

data = list(
  n = nrow(X),
  p = ncol(X),
  g = nrow(X),
  y = as.vector(y),
  X = X,
  mu_Sigma = Sigma,
  mu_b = rep(0, 4)
)

WIDTH = 3000
HEIGHT = 2400
RES = 300

# Does not includes model 7, which uses b_sd, fitted below.
for (i in c(1:6, 8)) {
  cat("Sampling and plotting model", i, "...\n")
  
  fit = sampling(
    readRDS(here::here("models", paste0("model_", i, ".rds"))), 
    data = data,
    refresh = 0
  )
  if ("sigma_p" %in% names(fit)) {
    # Prior
    png(
      here("tests/plots", paste0("prior_m", i, ".png")), 
      width = WIDTH, height = HEIGHT, res = RES
    )
    pairs(fit, pars = c("beta_p", "sigma_p"), main = paste0("Prior M", i))
    dev.off() 
  }
  # Posterior
  png(
    here("tests/plots", paste0("posterior_m", i, ".png")), 
    width = WIDTH, height = HEIGHT, res = RES
  )
  pairs(fit, pars = c("beta", "sigma"), main = paste0("Posterior M", i))
  dev.off() 
  cat("Done!\n")
}


# Model 7, b_sd = 1
data$b_sd = 1
fit = sampling(
  readRDS(here::here("models/model_7.rds")), 
  data = data,
  refresh = 0
)

# Prior
png(
  here("tests/plots/prior_m7_b_sd_1.png"), 
  width = WIDTH, height = HEIGHT, res = RES
)
pairs(fit, pars = c("beta_p", "sigma_p"), main = "Prior M7, b_sd = 1")
dev.off() 


# Posterior
png(
  here("tests/plots/posterior_m7_b_sd_1.png"), 
  width = WIDTH, height = HEIGHT, res = RES
)
pairs(fit, pars = c("beta", "sigma"), main = "Posterior M7, b_sd = 1")
dev.off() 

# Model 7, b_sd = 2.5
data$b_sd = 2.5
fit = sampling(
  readRDS(here::here("models/model_7.rds")), 
  data = data,
  refresh = 0
)

# Prior
png(
  here("tests/plots/prior_m7_b_sd_2_5.png"), 
  width = WIDTH, height = HEIGHT, res = RES
)
pairs(fit, pars = c("beta_p", "sigma_p"), main = "Prior M7, b_sd = 2.5")
dev.off() 


# Posterior
png(
  here("tests/plots/posterior_m7_b_sd_2_5.png"), 
  width = WIDTH, height = HEIGHT, res = RES
)
pairs(fit, pars = c("beta", "sigma"), main = "Posterior M7, b_sd = 2.5")
dev.off()
