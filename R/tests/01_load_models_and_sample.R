library(rstan)

# Toy sample --------------------------------------------------------------
SIZE = 20
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

# g = n, log(sigma) ~ uniform
fit = sampling(
  readRDS(here::here("models", "model_1.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 1")

# g = n, sigma ~ exp(1 / sd(y))
fit = sampling(
  readRDS(here::here("models", "model_2.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 2")

# g = n, sigma ~ sigma ~ student_t(nu = 4, sigma = sd(y))
fit = sampling(
  readRDS(here::here("models", "model_3.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 3")

# From here we don't have fixed 'g'.
data[["g"]] = NULL

# g ~ student_t(3, 0, 3), log(sigma) ~ uniform
fit = sampling(
  readRDS(here::here("models", "model_4.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 4")

# g ~ student_t(3, 0, 3), sigma ~ exp(1 / sd(y))
fit = sampling(
  readRDS(here::here("models", "model_5.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 5")

# g ~ student_t(3, 0, 3), sigma ~ student_t(nu = 4, sigma = sd(y))
fit = sampling(
  readRDS(here::here("models", "model_6.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 6")

# Independent N(0, 1) priors for all, sigma ~ exp(1 / sd(y))
fit = sampling(
  readRDS(here::here("models", "model_7.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 7")

# Independent N(0, 2.5) priors for all, sigma ~ exp(1 / sd(y))
fit = sampling(
  readRDS(here::here("models", "model_8.rds")), 
  data = data,
  refresh = 0
)
pairs(fit, pars = c("beta", "sigma"), main = "Model 8")
