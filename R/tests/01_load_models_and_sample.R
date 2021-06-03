library(rstan)

# Toy sample --------------------------------------------------------------
SIZE = 25
set.seed(1995)

x1 <- rnorm(SIZE)
x2 <- 0.35 + 0.15 * x1 + rnorm(SIZE)
x3 <- 0.5 + 0.3 * x1 + 0.2 * x2 + rnorm(SIZE)
x4 <- -0.7 + 0.2 * x1 + 0.1 * x3 + rnorm(SIZE)

X <- cbind(x1, x2, x3, x4)
X <- scale(X)
b_true <- c(2, 0.8, -1.5, -0.3)
sigma_true <- 2
y <- X %*% b_true + rnorm(SIZE, sd=sigma_true)
g <- nrow(X)
Sigma <- solve(t(X) %*% X)

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

for (i in 1:8) {
  cat("Sampling and plotting model", i, "...\n")
  fit = sampling(
    readRDS(here::here("models", paste0("model_", i, ".rds"))), 
    data = data,
    refresh = 0
  )
  # Prior
  png(
    here::here("tests", "plots", paste0("prior_m", i, ".png")), 
    width = WIDTH, height = HEIGHT, res = RES
  )
  pairs(fit, pars = c("beta_p", "sigma_p"), main = paste0("Prior M", i))
  dev.off() 
  # Posterior
  png(
    here::here("tests", "plots", paste0("posterior_m", i, ".png")), 
    width = WIDTH, height = HEIGHT, res = RES
  )
  pairs(fit, pars = c("beta", "sigma"), main = paste0("Posterior M", i))
  dev.off() 
  cat("Done!\n")
}
  