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
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | X[i] * beta, sigma);
  }
}
