data {
  
  // Metadata
  int N;             // Number of participants
  int T;             // Number of trials
  
  // Indices
  array[T] int subj_ix;    // Participant number for each datapoint
  
  // Data
  array[T] real Y;             // Dependent variable: self-reported mood valence
  array[T] int  trial_num;     // Animal-specific trial number (0-indexed)
  
}

parameters {
  
  // Group-level means
  real w_0_mu_pr;
  real response_sd_mu_pr;
  
  // Group-level SDs
  real<lower=0> w_0_sigma_pr;
  real<lower=0> response_sd_sigma_pr;
  
  // Participant-level parameters
  vector[N] w_0_pr;
  vector[N] response_sd_pr;
  
}

transformed parameters {
  
  vector[N] w_0;
  vector[N] response_sd;
  
  for (p_ix in 1:N){
    w_0[p_ix]              = w_0_mu_pr + w_0_sigma_pr * w_0_pr[p_ix];
    response_sd[p_ix]      = exp(response_sd_mu_pr + response_sd_sigma_pr * response_sd_pr[p_ix]);
  }
}

model {
  
  // group-level priors for means
  w_0_mu_pr              ~ normal(0, 1);
  response_sd_mu_pr      ~ normal(0, 10);
  
  // group-level priors for SDs
  w_0_sigma_pr              ~ exponential(0.1);
  response_sd_sigma_pr      ~ exponential(0.1);
  
  // subject-level priors
  w_0_pr              ~ normal(0, 1);
  response_sd_pr      ~ normal(0, 1);

  // predicted mean and sd (and their transformed cousins mu, phi, A, and B)
  vector[T] pred_mean = rep_vector(0, T);
  vector[T] mu;
  vector[T] phi;
  vector[T] A;
  vector[T] B;
  
  // loop over trials
  for (trial_ix in 1:T){
    
    // predicted mean - intercept
    pred_mean[trial_ix] += w_0[subj_ix[trial_ix]];
    
    // transform for beta regression
    mu[trial_ix]  = inv_logit(pred_mean[trial_ix]);
    phi[trial_ix] = response_sd[subj_ix[trial_ix]];
    A[trial_ix]   = mu[trial_ix] * phi[trial_ix] + machine_precision(); // add machine precision to prevent underflow (A = 0)
    B[trial_ix]   = phi[trial_ix] - (mu[trial_ix] * phi[trial_ix]) + machine_precision();  // add machine precision to prevent underflow (B = 0)
    
  }
  
  // responses distributed as beta
  Y ~ beta(A, B);
  
}

generated quantities {

  // predicted mean and sd (and their transformed cousins mu, phi, A, and B)
  vector[T] pred_mean = rep_vector(0, T);
  vector[T] mu;
  vector[T] phi;
  vector[T] A;
  vector[T] B;
  
  // containers for log likelihood and predicted choice
  array[T] real log_likelihood;
  array[T] real mood_pred;
  
  // loop over trials
  for (trial_ix in 1:T){
    
    // predicted mean - intercept
    pred_mean[trial_ix] += w_0[subj_ix[trial_ix]];
    
    // transform for beta regression
    mu[trial_ix]  = inv_logit(pred_mean[trial_ix]);
    phi[trial_ix] = response_sd[subj_ix[trial_ix]];
    A[trial_ix]   = mu[trial_ix] * phi[trial_ix] + machine_precision(); // add machine precision to prevent underflow (A = 0)
    B[trial_ix]   = phi[trial_ix] - (mu[trial_ix] * phi[trial_ix]) + machine_precision();  // add machine precision to prevent underflow (B = 0)
    
    // log-likelihood of observations
    log_likelihood[trial_ix] = beta_lpdf(Y[trial_ix] | A[trial_ix], B[trial_ix]);
    
  }
  
  // predicted mood ratings
  mood_pred = beta_rng(A, B);
  
}
