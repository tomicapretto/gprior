
* Model 1: g = n and log(sigma) ~ uniform
* Model 2: g = n and sigma ~ exponential(1/sd(y))
* Model 3: g = n and sigma ~ student_t(nu=4, sigma=sd(y))
* Model 4: g ~ student_t(3, 0, 3) and log(sigma) ~ uniform
* Model 5: g ~ student_t(3, 0, 3) and sigma ~ exponential(1/sd(y))
* Model 6: g ~ student_t(3, 0, 3) and student_t(nu=4, sigma=sd(y))
* Model 7: No g. Betas have independent N(0, 1) and sigma ~ exponential(1/sd(y))
* Model 8: No g. Betas have independent N(0, 2.5) and sigma ~ exponential(1/sd(y))


Notes

* Predictors are standardized
* student_t(3, 0, 3) was chosen so P_95 ~ 10, which means we inflate b_ols covariance matrix by at most 10 for the prior.
* sigma ~ student_t(nu=4, sigma=sd(y)) is what is currently done in Bambi.
* sigma ~ exponential(1/sd(y)) currently done in rstanarm.