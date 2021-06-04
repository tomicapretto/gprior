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
  y ~ normal(to_vector(X * beta), sigma);
}
