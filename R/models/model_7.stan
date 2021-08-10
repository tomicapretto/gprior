data {
  int n;
  int p;
  real<lower=0> b_sd;
  vector[n] y;
  vector[p] mu_b;
  matrix[n, p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;
}
model {
  beta ~ normal(mu_b, b_sd * sd(y));
  y ~ normal(to_vector(X * beta), sigma);
  sigma ~ exponential(1 / sd(y));
}
generated quantities {
  real<lower=0> sigma_p = exponential_rng(1 / sd(y));
  array[p] real beta_p = normal_rng(mu_b, b_sd * sd(y));
  array[n] real y_p = normal_rng(to_vector(X * to_vector(beta_p)), sigma_p);
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | X[i] * beta, sigma);
  }
}
