data {
  int<lower=0> N;  // Number of data points
  int<lower=0> J;  // Number of groups
  vector[N] x;     // Predictor variable
  vector[N] y;     // Outcome variable
  int<lower=1, upper=J> group[N]; // Group identifier for each data point
}

parameters {
  real alpha;                    // Global intercept
  real mu;                       // Mean slope across groups
  real<lower=0> tau;             // Standard deviation of slopes across groups
  vector[J] beta_raw;            // "Raw" slope for each group
  real<lower=0> sigma;           // Residual standard deviation
}

transformed parameters {
  vector[J] beta = mu + tau * beta_raw;  // "True" slope for each group
}

model {
  beta_raw ~ normal(0, 1);  // Prior for the non-centered parameter
  y ~ normal(alpha + beta[group] .* x, sigma);  // Likelihood
}
