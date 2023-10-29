data {
  int<lower=0> N;  // Number of data points
  int<lower=0> J;  // Number of groups
  array[N] real x;       // Predictor variable
  array[N] real y;       // Outcome variable
  array[N] int<lower=1, upper=J> group; // Group identifier for each data point
}

parameters {
  real alpha;              // Global intercept
  real mu;                 // Mean slope across groups
  real<lower=0> tau;       // Standard deviation of slopes across groups
  array[J] real beta;            // Slope for each group
  real<lower=0> sigma;     // Residual standard deviation
}

model {
  beta ~ normal(mu, tau);  // Group-level slopes
  
  for (i in 1:N) {
    y[i] ~ normal(alpha + beta[group[i]] * x[i], sigma);  // Likelihood
  }
}
