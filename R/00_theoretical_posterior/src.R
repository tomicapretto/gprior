# This script checks that the implementation of the classic g-prior in Stan
# matches the theoretical posterior.

# Libraries ---------------------------------------------------------------
library(ggplot2)
library(rstan)

# Helpers definition ------------------------------------------------------
get_beta_ols <- function(y, X) {
  solve((t(X) %*% X)) %*% (t(X) %*% y)
}

get_s_sq <- function(y, X, b_ols) {
  m <- y - X %*% b_ols
  t(m) %*% m
}

get_variance_posterior <- function(y, X, b_prior, b_ols, g, size = 1000) {
  s_sq <- get_s_sq(y, X, b_ols)
  b_diff <- b_prior - b_ols
  a <- nrow(X) / 2
  b <- 0.5 * (s_sq + (t(b_diff) %*% (t(X) %*% X) %*% b_diff) / (g + 1))
  return(1 / rgamma(size, a, b))
}

get_beta_posterior <- function(x, b_prior, b_ols, variance, g) {
  q <- g / (1 + g)
  mu <- q * b_ols + (1 - q) * b_prior
  Sigma <- q * solve(t(X) %*% X)
  l <- lapply(variance, function(var) mvtnorm::rmvnorm(1, mu, var * Sigma))
  values <- as.vector(do.call("rbind", l))
  data.frame(
    value = values,
    coef = rep(paste0("β", 1:length(mu)), each = length(values) / length(mu))
  )
}

# Setup -------------------------------------------------------------------
SIZE = 20
set.seed(1995)

x1 <- rnorm(SIZE)
x2 <- 0.35 + 0.15 * x1 + rnorm(SIZE)
x3 <- 0.5 + 0.3 * x1 + 0.2 * x2 + rnorm(SIZE)
x4 <- -0.7 + 0.2 * x1 + 0.1 * x3 + rnorm(SIZE)

X <- cbind(x1, x2, x3, x4)
X <- X - colMeans(X)
b_true <- c(2, 0.5, -1.5, -0.3)
sigma_true <- 2
y <- X %*% b_true + rnorm(SIZE, sd=sigma_true)
g <- nrow(X)
Sigma <- solve(t(X) %*% X)

# Stan model --------------------------------------------------------------
# Define model
model_str <- "data {
  int n; 
  int p; 
  real g;
  vector[n] y; 
  matrix[n, p] X;
  matrix[p, p] Sigma; // (X'X)^(-1)
  vector[p] mu;
}

parameters {
  vector[p] beta;
  real log_sigma_sq;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = sqrt(exp(log_sigma_sq));
}

model {
  log_sigma_sq ~ uniform(-5, 5);
  beta ~ multi_normal(mu, g * pow(sigma, 2) * Sigma);
  y ~ normal(to_vector(X * beta), sigma);
}
"

# Compile and sample
model <- stan_model(model_code = model_str)
model_fit <- sampling(
  model, 
  data = list(
    n = nrow(X),
    p = ncol(X),
    g = nrow(X),
    y = as.vector(y),
    X = X,
    Sigma = Sigma,
    mu = rep(0, 4)
  )
)

# Extract samples from the posterior
posterior <- extract(model_fit, c("beta", "sigma"))

# Convert to data.frame for usage with ggplot2
beta_posterior_trace <- data.frame(
  value = as.vector(posterior$beta),
  coef = rep(paste0("β", seq(ncol(posterior$beta))), each = nrow(posterior$beta)),
  type = "Stan"
)

sigma_posterior_trace <- data.frame(
  value = as.vector(posterior$sigma),
  coef = "σ",
  type = "Stan"
)

posterior_trace <- rbind(beta_posterior_trace, sigma_posterior_trace)

# Samples from theoretical posterior --------------------------------------
b_ols <- get_beta_ols(y, X)
sigma_posterior_exact <- data.frame(
  value = get_variance_posterior(y, X, rep(0, 4), b_ols, g) ^ 0.5,
  coef = "σ",
  type = "Exact"
)

beta_posterior_exact <- get_beta_posterior(X, rep(0, 4), b_ols, sigma_posterior_exact$value ^ 2, g)
beta_posterior_exact$type <- "Exact"

posterior_exact <- rbind(beta_posterior_exact, sigma_posterior_exact)

# Plots -------------------------------------------------------------------
colors <- c("#003f5c", "#7a5195", "#ef5675", "#ffa600")
ggplot(rbind(posterior_exact,posterior_trace)) + 
  geom_histogram(
    aes(x = value, y = ..density.., fill = type),
    position = "identity",
    bins = 50, 
    alpha = 0.5
  ) + 
  geom_point(
    aes(x = val), 
    y = 0,
    size = 2,
    color = "gray30",
    data = data.frame(
      val = c(b_true, sigma_true), 
      coef = c(paste0("β", 1:4), "σ")
    )
  ) + 
  labs(
    x = "Value",
    caption = "Point indicates true value"
  ) + 
  scale_fill_manual(values = colors[c(1, 4)]) + 
  facet_wrap(
    vars(coef), 
    nrow = 5, 
    scales = "free_x"
  ) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave("00_theoretical_posterior/posterior.png", height = 8, width = 5)
