data {
  int<lower=0> N;  // Number of participants
  int<lower=0> T;  // Total number of trials across all participants
  int<lower=0> S;  // Total number of sessions across all participants
  array[T] int<lower=1, upper=N> subj_ix; // Subject index for each trial
  array[T] int<lower=1, upper=S> session; // Session index for each trial
  array[T] int<lower=1> Y; // Ordinal outcome variable for each trial
  array[S] int<lower=1, upper=N> session_subj; // Subject index for each session
  vector[T] RPE;  // RPE for each trial
}
parameters {
  real w_0_mu_pr;  // Group-level mean prior for intercept
  real<lower=0> w_0_sigma_pr;  // Group-level standard deviation prior for intercept
  real b_session_mu_pr;  // Group-level mean prior for session effect
  real<lower=0> b_session_sigma_pr;  // Group-level standard deviation prior for session effect
  vector[N] w_0_pr;  // Individual-level random intercepts
  vector[S] b_session_pr;  // Session-level random intercepts
  ordered[4] cutpoints;  // Cutpoints for ordinal regression
  real beta_RPE;  // Weight for RPE
  real<lower=0, upper=1> gamma;  // Decay parameter for the exponential decay
}

transformed parameters {
  vector[N] w_0 = w_0_mu_pr + w_0_sigma_pr * w_0_pr;  // Individual-level intercepts
  vector[S] b_session = b_session_mu_pr + b_session_sigma_pr * b_session_pr;  // Session-level intercepts
  vector[T] weighted_RPE;  // Declare weighted_RPE here
  
  for (t in 1:T) {
    int subj = subj_ix[t];
    int sess = session[t];
    weighted_RPE[t] = 0;
    
    for (past_t in 1:t) {
      if (subj_ix[past_t] == subj && session[past_t] <= sess) {
        weighted_RPE[t] += pow(gamma, sess - session[past_t]) * RPE[past_t];
      }
    }
  }
}

model {
  w_0_mu_pr ~ normal(0, 1);
  w_0_sigma_pr ~ exponential(0.1);
  b_session_mu_pr ~ normal(0, 1);
  b_session_sigma_pr ~ exponential(0.1);
  beta_RPE ~ normal(0, 1);  // Prior for RPE
  gamma ~ beta(1, 1);  // Prior for gamma, uniform between 0 and 1
  w_0_pr ~ normal(0, 1);
  b_session_pr ~ normal(0, 1);
  cutpoints ~ normal(0, 1);

  vector[T] pred_mean;

  for (t in 1:T) {
    int subj = subj_ix[t];
    int sess = session[t];

    // Removed the nested loop for calculating weighted_RPE as it's already calculated in the transformed parameters block

    pred_mean[t] = w_0[subj] + b_session[sess] + beta_RPE * weighted_RPE[t];
  }

  Y ~ ordered_logistic(pred_mean, cutpoints);
}

generated quantities {
  array[T] int Y_pred;  // Predicted categories for Y
  vector[T] log_lik;  // Pointwise log-likelihood
  vector[T] pred_mean = w_0[subj_ix] + b_session[session] + beta_RPE * weighted_RPE;
  
  for (t in 1:T) {
    log_lik[t] = ordered_logistic_lpmf(Y[t] | pred_mean[t], cutpoints);
    Y_pred[t] = ordered_logistic_rng(pred_mean[t], cutpoints);  // Generating predicted categories for Y
  }
}

