# Do we scale all predictors?
# Do we include intercept?


# Remember
# real<lower=0> sigma;
# real<lower=0> g;

# Need to convert mu and sigma to alpha and beta in gamma.

# All combinations between
# g = n
# g ~ gamma("g", mu=n ^ 0.8, sigma=5)
# g ~ student_t(3, 0, 1)


# AND

# log(sigma) ~ uniform
# sigma ~ exp(1) --> check if this is after scaling "y".
# sigma ~ student_t(nu=4, sigma=sd(y))


# And finally,
# Independent N(0, 1) priors for all, except intercepts who have N(mean, 1).
# Aki prefers student_t(3,0,1),



# Things I can't miss -----------------------------------------------------
# 1. Evaluate prior in terms of the distribution on the response scale.
#    What do we want? Uniform? Bell-shaped? Centered around which value?
# 2. Evaluate coverage of 95% probability intervals.
# 3. Measure sampling time. 
# 4. Count divergences.
# 5. Effective sample size.
# 6. Store R-hat.
# 7. Models with low effects, medium effects, and large effects in terms of SD.
#    Something like < 0.5, > 0.5 & < 2.5, and > 2.5 (with some >> 2.5)
# 8. Implied mean standard deviation for priors of the coefficients.



# Here's an idea for not getting tripped up with default priors: For each 
# parameter (or other qoi), compare the posterior sd to the prior sd. 
# If the posterior sd for any parameter (or qoi) is more than 0.1 times the 
# prior sd, then print out a note: "The prior distribution for this parameter 
# is informative." Then the user can go back and check that the default prior 
# makes sense for this particular example.


# One principle: write down what you think the prior should be, then spread it out.
# The idea is that the cost of setting the prior too narrow is more severe than 
# the cost of setting it too wide. I've been having trouble formalizing this idea.


# One principle: write down what you think the prior should be, then spread it out. 
# The idea is that the cost of setting the prior too narrow is more severe than 
# the cost of setting it too wide. I've been having trouble formalizing this idea.