data {
  int<lower=0> N;  // Number of observations
  int<lower=1> T;  // Number of trials
  vector[N] RPE; 
  vector[N] outcome; 
  vector[N] stimulus;
  vector[N] zmoodpre;
  vector[N] zcontrol;
  vector[N] trial;
  vector[N] happiness;
}

parameters {
  real alpha;  // Intercept
  vector[6] beta;  // Slopes
  real<lower=0, upper=1> gamma;  // Decay parameter
}

model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  gamma ~ beta(1, 1);  // Prior for gamma

  vector[N] mu_happiness;
  for (n in 1:N) {
    // Compute weighted values with exponential decay
    real weighted_outcome_local = pow(gamma, trial[n]-1) * outcome[n];
    real weighted_stimulus_local = pow(gamma, trial[n]-1) * stimulus[n];
    real weighted_RPE_local = pow(gamma, trial[n]-1) * RPE[n];

    // Compute the mean of the Gaussian for each observation
    mu_happiness[n] = alpha + 
      beta[1] * weighted_outcome_local +
      beta[2] * weighted_stimulus_local +
      beta[3] * weighted_RPE_local +
      beta[4] * zmoodpre[n] +
      beta[5] * zcontrol[n] +
      beta[6] * trial[n];
  }
  happiness ~ normal(mu_happiness, 1);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for LOO
  vector[N] y_pred;   // Predicted values

  for (n in 1:N) {
    // Compute weighted values with exponential decay
    real weighted_outcome_local = pow(gamma, trial[n]-1) * outcome[n];
    real weighted_stimulus_local = pow(gamma, trial[n]-1) * stimulus[n];
    real weighted_RPE_local = pow(gamma, trial[n]-1) * RPE[n];

    // Compute the mean of the Gaussian for each observation
    real mu_happiness = alpha + 
      beta[1] * weighted_outcome_local +
      beta[2] * weighted_stimulus_local +
      beta[3] * weighted_RPE_local +
      beta[4] * zmoodpre[n] +
      beta[5] * zcontrol[n] +
      beta[6] * trial[n];

    log_lik[n] = normal_lpdf(happiness[n] | mu_happiness, 1);
    y_pred[n] = normal_rng(mu_happiness, 1);
  }
}
