data {
  int<lower=0> T;  // Total number of trials
  int<lower=0> N;  // Number of participants
  array[T] int<lower=1, upper=N> subj_ix; // Subject index for each trial
  vector[T] Y;  // Standardized happiness for each trial
  vector[T] RPE;  // RPE for each trial
  vector[T] outcome; // Outcome for each trial
  vector[T] stimulus; // Stimulus for each trial
  vector[T] zmoodpre; // zmoodpre for each trial
  vector[T] zcontrol; // zcontrol for each trial
  vector[T] trial; // Trial for each trial
}

parameters {
  real w_0;  // intercept
  real beta_RPE;  // Weight for RPE
  real beta_outcome; // Weight for outcome
  real beta_stimulus; // Weight for stimulus
  real beta_zmoodpre; // Weight for zmoodpre
  real beta_zcontrol; // Weight for zcontrol
  real beta_trial; // Weight for trial
  real<lower=0, upper=1> gamma;  // Decay parameter for the exponential decay
  real<lower=0> sigma;  // Residual standard deviation
}

transformed parameters {
  vector[T] weighted_RPE = rep_vector(0, T);
  vector[T] weighted_outcome = rep_vector(0, T);
  vector[T] weighted_stimulus = rep_vector(0, T);
  
  for (t in 1:T) {
    for (past_t in 1:t) {
      weighted_RPE[t] += pow(gamma, (t - past_t)) * RPE[past_t];
      weighted_outcome[t] += pow(gamma, (t - past_t)) * outcome[past_t];
      weighted_stimulus[t] += pow(gamma, (t - past_t)) * stimulus[past_t];
    }
  }
}

model {
  w_0 ~ normal(0, 1);
  beta_RPE ~ normal(0, 1);
  beta_outcome ~ normal(0, 1);
  beta_stimulus ~ normal(0, 1);
  beta_zmoodpre ~ normal(0, 1);
  beta_zcontrol ~ normal(0, 1);
  beta_trial ~ normal(0, 1);
  gamma ~ beta(1, 1);
  sigma ~ exponential(0.1);
  
  Y ~ normal(w_0 + 
             beta_RPE * weighted_RPE +
             beta_outcome * weighted_outcome +
             beta_stimulus * weighted_stimulus +
             beta_zmoodpre * zmoodpre +
             beta_zcontrol * zcontrol +
             beta_trial * trial, sigma);
}

generated quantities {
  vector[T] Y_pred;  // Predicted values for Y
  vector[T] log_lik;  // Pointwise log-likelihood
  
  for (t in 1:T) {
    // Compute the predicted mean for each observation
    Y_pred[t] = normal_rng(
                  w_0 + 
                  beta_RPE * weighted_RPE[t] +
                  beta_outcome * weighted_outcome[t] +
                  beta_stimulus * weighted_stimulus[t] +
                  beta_zmoodpre * zmoodpre[t] +
                  beta_zcontrol * zcontrol[t] +
                  beta_trial * trial[t], sigma);
    
    // Compute the pointwise log-likelihood
    log_lik[t] = normal_lpdf(Y[t] | Y_pred[t], sigma);
  }
}

