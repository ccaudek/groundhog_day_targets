//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=0> N; // total number of observations
  int<lower=0> J; // total number of participants
  array[N] int<lower=1, upper=J> participant; // participant index for each observation
  
  array[N] real outcome;
  array[N] real stimulus;
  array[N] real RPE;
  array[N] real zmoodpre;
  array[N] real zcontrol;
  array[N] real trial;
  array[N] real happiness;
}

parameters {
  real<lower=0> sigma; // residual error
  
  // Non-centered participant-level parameters
  array[J] real w0_raw;
  array[J] real w1_raw;
  array[J] real w2_raw;
  array[J] real w3_raw;
  array[J] real w4_raw;
  array[J] real w5_raw;
  array[J] real w6_raw;
  array[J] real gamma_raw;

  // Hyperparameters
  real mu_w0;
  real mu_w1;
  real mu_w2;
  real mu_w3;
  real mu_w4;
  real mu_w5;
  real mu_w6;
  real mu_gamma;

  real<lower=0> sigma_w0;
  real<lower=0> sigma_w1;
  real<lower=0> sigma_w2;
  real<lower=0> sigma_w3;
  real<lower=0> sigma_w4;
  real<lower=0> sigma_w5;
  real<lower=0> sigma_w6;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  // Participant-level parameters, transformed
  array[J] real w0;
  array[J] real w1;
  array[J] real w2;
  array[J] real w3;
  array[J] real w4;
  array[J] real w5;
  array[J] real w6;
  array[J] real gamma;

  for (j in 1:J) {
    w0[j] = mu_w0 + sigma_w0 * w0_raw[j];
    w1[j] = mu_w1 + sigma_w1 * w1_raw[j];
    w2[j] = mu_w2 + sigma_w2 * w2_raw[j];
    w3[j] = mu_w3 + sigma_w3 * w3_raw[j];
    w4[j] = mu_w4 + sigma_w4 * w4_raw[j];
    w5[j] = mu_w5 + sigma_w5 * w5_raw[j];
    w6[j] = mu_w6 + sigma_w6 * w6_raw[j];
    gamma[j] = mu_gamma + sigma_gamma * gamma_raw[j];
  }
}


model {
  // Priors for residual error
  sigma ~ normal(0, 1);

  // Priors for non-centered raw parameters
  w0_raw ~ std_normal();
  w1_raw ~ std_normal();
  w2_raw ~ std_normal();
  w3_raw ~ std_normal();
  w4_raw ~ std_normal();
  w5_raw ~ std_normal();
  w6_raw ~ std_normal();
  gamma_raw ~ std_normal();
  
  // Priors for hyperparameters
  mu_w0 ~ normal(0, 1);
  mu_w1 ~ normal(0, 1);
  mu_w2 ~ normal(0, 1);
  mu_w3 ~ normal(0, 1);
  mu_w4 ~ normal(0, 1);
  mu_w5 ~ normal(0, 1);
  mu_w6 ~ normal(0, 1);
  mu_gamma ~ normal(0, 1);

  sigma_w0 ~ normal(0, 1);
  sigma_w1 ~ normal(0, 1);
  sigma_w2 ~ normal(0, 1);
  sigma_w3 ~ normal(0, 1);
  sigma_w4 ~ normal(0, 1);
  sigma_w5 ~ normal(0, 1);
  sigma_w6 ~ normal(0, 1);
  sigma_gamma ~ normal(0, 1);

  // Likelihood
  for (i in 1:N) {
    int idx = participant[i];
    real weighted_outcome = 0;
    real weighted_stimulus = 0;
    real weighted_RPE = 0;
    
    for (t in 1:i) {
      weighted_outcome += pow(gamma[idx], i - t) * outcome[t];
      weighted_stimulus += pow(gamma[idx], i - t) * stimulus[t];
      weighted_RPE += pow(gamma[idx], i - t) * RPE[t];
    }
    
    happiness[i] ~ normal(
      w0[idx] + 
      w1[idx] * weighted_outcome + 
      w2[idx] * weighted_stimulus + 
      w3[idx] * weighted_RPE + 
      w4[idx] * zmoodpre[i] + 
      w5[idx] * zcontrol[i] + 
      w6[idx] * trial[i], 
      sigma
    );
  }
}



