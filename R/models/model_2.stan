data {
  int n;
  int p;
  real g;
  vector[n] y;
  vector[p] mu_b;
  matrix[p, p] mu_Sigma;
  matrix[n, p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;
}
model {
  beta ~ multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma);
  y ~ normal(to_vector(X * beta), sigma);
  sigma ~ exponential(1 / sd(y));
}
generated quantities {
  real sigma_p = exponential_rng(1 / sd(y));
  vector[p] beta_p = multi_normal_rng(mu_b, g * pow(sigma_p, 2) * mu_Sigma);
  real y_p[n] = normal_rng(to_vector(X * beta_p), sigma_p);
}
