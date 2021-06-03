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
  real log_sigma_sq;
  real<lower=0> g;
}
transformed parameters {
  real<lower=0> sigma = sqrt(exp(log_sigma_sq));
}
model {
  beta ~ multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma);
  y ~ normal(to_vector(X * beta), sigma);
  log_sigma_sq ~ uniform(-5, 5);
  g ~ student_t(3, 0, 4);
}
generated quantities {
  real log_sigma_sq_p = uniform_rng(-5, 5);
  real<lower=0> sigma_p = sqrt(exp(log_sigma_sq_p));
  real g_p = fabs(student_t_rng(3, 0, 4));
  vector[p] beta_p = multi_normal_rng(mu_b, g_p * pow(sigma_p, 2) * mu_Sigma);
  real y_p[n] = normal_rng(to_vector(X * beta_p), sigma_p);
}
