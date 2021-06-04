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
  g ~ student_t(3, 0, 4);
}
