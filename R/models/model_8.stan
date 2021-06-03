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
