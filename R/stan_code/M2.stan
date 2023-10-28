//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//

data {
  int N; // Number of participants
  int T; // Number of trials
  
  array[T] int subj_ix; // Participant index for each datapoint
  
  array[T] real happiness; // Dependent variable: momentary happiness
  array[T] real outcome;
  array[T] real stimulus;
  array[T] real RPE;
  array[T] real zmoodpre;
  array[T] real zcontrol;
  array[T] int trial_num;
}

parameters {
  // Group-level means
  real w_0_mu;
  real gamma_mu;
  
  // Group-level SDs
  real<lower=0> w_0_sigma;
  real<lower=0> gamma_sigma;
  
  // Participant-level parameters (non-centered)
  vector[N] w_0_nc;
  vector[N] gamma_nc;
}

transformed parameters {
  vector[N] w_0;
  vector[N] gamma;
  
  for (p in 1:N){
    w_0[p] = w_0_mu + w_0_sigma * w_0_nc[p];
    gamma[p] = inv_logit(gamma_mu + gamma_sigma * gamma_nc[p]);
  }
}

model {
  // Group-level priors
  w_0_mu ~ normal(0, 1);
  gamma_mu ~ normal(0, 1);
  
  w_0_sigma ~ exponential(0.1);
  gamma_sigma ~ exponential(0.1);
  
  // Non-centered individual-level priors
  w_0_nc ~ normal(0, 1);
  gamma_nc ~ normal(0, 1);
  
  for (trial in 1:T) {
    int idx = subj_ix[trial];
    real weighted_outcome = 0;
    real weighted_stimulus = 0;
    real weighted_RPE = 0;
    
    for (t in 1:trial) {
      weighted_outcome += pow(gamma[idx], trial - t) * outcome[t];
      weighted_stimulus += pow(gamma[idx], trial - t) * stimulus[t];
      weighted_RPE += pow(gamma[idx], trial - t) * RPE[t];
    }
    
    real mu = fmax(0.001, inv_logit(
      w_0[idx] + 
      weighted_outcome + 
      weighted_stimulus + 
      weighted_RPE + 
      zmoodpre[trial] + 
      zcontrol[trial] + 
      trial_num[trial]
    ));
    real phi = fmax(0.001, 1); // Dispersion parameter, you can model this too if needed
    
    happiness[trial] ~ beta(mu * phi, (1 - mu) * phi);
  }
}

generated quantities {
  array[T] real happiness_pred;  // Predicted momentary happiness
  array[T] real log_lik;         // Log-likelihood for each observation
  
  for (trial in 1:T) {
    int idx = subj_ix[trial];
    real weighted_outcome = 0;
    real weighted_stimulus = 0;
    real weighted_RPE = 0;
    
    for (t in 1:trial) {
      weighted_outcome += pow(gamma[idx], trial - t) * outcome[t];
      weighted_stimulus += pow(gamma[idx], trial - t) * stimulus[t];
      weighted_RPE += pow(gamma[idx], trial - t) * RPE[t];
    }
    
    real mu = inv_logit(
      w_0[idx] + 
      weighted_outcome + 
      weighted_stimulus + 
      weighted_RPE + 
      zmoodpre[trial] + 
      zcontrol[trial] + 
      trial_num[trial]
    );
    real phi = 1; // Dispersion parameter, you can model this too if needed
    
    // Generate predicted momentary happiness
    happiness_pred[trial] = beta_rng(mu * phi, (1 - mu) * phi);
    
    // Compute log-likelihood for each observed value
    log_lik[trial] = beta_lpdf(happiness[trial] | mu * phi, (1 - mu) * phi);
  }
}


