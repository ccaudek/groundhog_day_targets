data {
  int<lower=0> N;  // Number of data points
  vector[N] alpha; // Intercept from individual models
  vector[N] beta;  // Coefficient for feedback from individual models
}

parameters {
  real mu_alpha;      // Population mean for alpha
  real<lower=0> sigma_alpha; // Population standard deviation for alpha
  real mu_beta;       // Population mean for beta
  real<lower=0> sigma_beta;  // Population standard deviation for beta
}

model {
  // Priors
  mu_alpha ~ normal(0, 10);
  sigma_alpha ~ normal(0, 10);
  mu_beta ~ normal(0, 10);
  sigma_beta ~ normal(0, 10);

  // Likelihood
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
}
