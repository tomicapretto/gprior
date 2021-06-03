data {
  int n;
  int p;
  vector[n] y;
  vector[p] mu_b;
  matrix[n, p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;
}
model {
  beta ~ normal(mu_b, 2.5);
  y ~ normal(to_vector(X * beta), sigma);
  sigma ~ exponential(1 / sd(y));
}
generated quantities {
  real<lower=0> sigma_p = exponential_rng(1/sd(y));
  real beta_p[p] = normal_rng(mu_b, 2.5);
  real y_p[n] = normal_rng(to_vector(X * to_vector(beta_p)), sigma_p);
}
