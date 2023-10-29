data {
  int<lower=0> N;  // Number of data points
  vector[N] x;     // Predictor variable
  vector[N] y;     // Outcome variable
}

parameters {
  real alpha;                 // Intercept
  real beta_raw;              // "Raw" slope parameter (non-centered)
  real<lower=0> sigma;        // Standard deviation
}

transformed parameters {
  real beta = 0 + 1 * beta_raw;  // "True" slope, shifted and scaled
}

model {
  beta_raw ~ normal(0, 1);  // Prior for the non-centered parameter
  y ~ normal(alpha + beta * x, sigma);
}
