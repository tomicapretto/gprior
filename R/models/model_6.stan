data {
  int n;
  int p;
  vector[n] y;
  vector[p] mu_b;
  matrix[p, p] mu_Sigma;
  matrix[n, p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> g;
  real<lower=0> sigma;
}
model {
  beta ~ multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma);
  y ~ normal(to_vector(X * beta), sigma);
  g ~ student_t(3, 0, 4);
  sigma ~ student_t(4, 0, sd(y));
}
generated quantities {
  real<lower=0> sigma_p = fabs(student_t_rng(4, 0, sd(y)));
  real<lower=0> g_p = fabs(student_t_rng(3, 0, 4));
  vector[p] beta_p = multi_normal_rng(mu_b, g_p * pow(sigma_p, 2) * mu_Sigma);;
  array[n] real y_p = normal_rng(to_vector(X * beta_p), sigma_p);
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | X[i] * beta, sigma);
  }
}
