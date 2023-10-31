data {
  int<lower=0> N; // Number of participants
  int<lower=0> T; // Total number of trials across all participants
  int<lower=0> S; // Total number of sessions across all participants
  array[T] int<lower=1, upper=N> subj_ix; // Subject index for each trial
  array[T] int<lower=1, upper=S> session; // Session index for each trial
  array[T] int<lower=1> Y; // Ordinal outcome variable for each trial
  array[S] int<lower=1, upper=N> session_subj; // Subject index for each session
  vector[T] RPE; // RPE for each trial
}
parameters {
  real w_0_mu_pr; // Group-level mean prior for intercept
  real<lower=0> w_0_sigma_pr; // Group-level standard deviation prior for intercept
  real b_session_mu_pr; // Group-level mean prior for session effect
  real<lower=0> b_session_sigma_pr; // Group-level standard deviation prior for session effect
  vector[N] w_0_pr; // Individual-level random intercepts
  vector[S] b_session_pr; // Session-level random intercepts
  ordered[4] cutpoints; // Cutpoints for ordinal regression
  real beta_RPE; // Weight for RPE
  real<lower=0, upper=1> gamma; // Decay parameter for the exponential decay
}
transformed parameters {
  vector[N] w_0; // Individual-level intercepts
  vector[S] b_session; // Session-level intercepts
  
  // Transforming the participant-level random intercepts
  for (n in 1 : N) {
    w_0[n] = w_0_mu_pr + w_0_sigma_pr * w_0_pr[n];
  }
  
  // Transforming the session-level random intercepts
  for (s in 1 : S) {
    b_session[s] = b_session_mu_pr + b_session_sigma_pr * b_session_pr[s];
  }
}
model {
  // Group-level priors
  w_0_mu_pr ~ normal(0, 1);
  w_0_sigma_pr ~ exponential(0.1);
  b_session_mu_pr ~ normal(0, 1);
  b_session_sigma_pr ~ exponential(0.1);
  beta_RPE ~ normal(0, 1); // Prior for RPE
  gamma ~ beta(1, 1); // Prior for gamma, uniform between 0 and 1
  
  // Participant-level priors
  w_0_pr ~ normal(0, 1);
  
  // Session-level priors
  b_session_pr ~ normal(0, 1);
  
  // Cutpoint priors
  cutpoints ~ normal(0, 1);
  
  // Ordered logistic likelihood with session and RPE from the current trial
  real weighted_RPE; // Declare the weighted_RPE variable outside the loop
  
  for (t in 1 : T) {
    weighted_RPE = 0;
    int subj = subj_ix[t];
    int sess = session[t];
    
    // Calculate the exponentially decayed sum of past RPEs for the current subject
    for (past_t in 1 : t) {
      if (subj_ix[past_t] == subj && session[past_t] <= sess) {
        weighted_RPE += pow(gamma, sess - session[past_t]) * RPE[past_t];
      }
    }
    
    // Predicted mean includes the weighted RPE and random effect for session
    real pred_mean = w_0[subj] + b_session[sess] + beta_RPE * weighted_RPE;
    
    // Ordered logistic model
    Y[t] ~ ordered_logistic(pred_mean, cutpoints);
  }
}

generated quantities {
  vector[T] log_lik;  // Pointwise log-likelihood
  
  for (t in 1:T) {
    real weighted_RPE = 0;
    int subj = subj_ix[t];
    int sess = session[t];
    
    // Calculate the exponentially decayed sum of past RPEs for the current subject
    for (past_t in 1:t) {
      if (subj_ix[past_t] == subj && session[past_t] <= sess) {
        weighted_RPE += pow(gamma, sess - session[past_t]) * RPE[past_t];
      }
    }
    
    // Predicted mean includes the weighted RPE and random effect for session
    real pred_mean = w_0[subj] + b_session[sess] + beta_RPE * weighted_RPE;
    
    // Store pointwise log-likelihood for later analysis
    log_lik[t] = ordered_logistic_lpmf(Y[t] | pred_mean, cutpoints);
  }
}

