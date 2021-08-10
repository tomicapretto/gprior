# for the data, cv prediction error with loo::loo()
library(here)
library(loo)
library(rstan)

i=2

SIZE = 1000
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

fit = sampling(
  readRDS(here::here("models", paste0("model_", 7, ".rds"))), 
  data = append(data, list(b_sd = 2.5)),
  refresh = 0
)

# Extract pointwise log-likelihood
# using merge_chains=FALSE returns an array, which is easier to 
# use with relative_eff()
log_lik_1 <- extract_log_lik(fit, merge_chains = FALSE)

# as of loo v2.0.0 we can optionally provide relative effective sample sizes
# when calling loo, which allows for better estimates of the PSIS effective
# sample sizes and Monte Carlo error
r_eff <- relative_eff(exp(log_lik_1), cores = 2)

# preferably use more than 2 cores (as many cores as possible)
# will use value of 'mc.cores' option if cores is not specified
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
