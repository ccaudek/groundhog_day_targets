// Hierarchical fixed-effect model for momentary happiness.
//
// Sun Oct 29 10:57:12 CET 2023 
data {
  int<lower=0> N;  // Number of participants
  int<lower=0> T;  // Total number of trials across all participants
  int<lower=0> S;  // Total number of sessions across all participants
  array[T] int<lower=1, upper=N> subj_ix;  // Subject index for each trial
  array[T] int<lower=1, upper=S> session;  // Session index for each trial
  array[T] int<lower=1> Y;  // Ordinal outcome variable for each trial
  array[S] int<lower=1, upper=N> session_subj;  // Subject index for each session
}

parameters {
  real w_0_mu_pr;  // Group-level mean prior for intercept
  real<lower=0> w_0_sigma_pr;  // Group-level standard deviation prior for intercept
  real beta_session;  // Fixed effect for session
  vector[N] w_0_pr;  // Individual-level random intercepts
  ordered[4] cutpoints;  // Cutpoints for ordinal regression
}

transformed parameters {
  vector[N] w_0;  // Individual-level intercepts

  // Transforming the participant-level random intercepts
  for (n in 1:N) {
    w_0[n] = w_0_mu_pr + w_0_sigma_pr * w_0_pr[n];
  }
}

model {
  // Group-level priors
  w_0_mu_pr ~ normal(0, 1);
  w_0_sigma_pr ~ exponential(0.1);

  // Fixed effect prior
  beta_session ~ normal(0, 1);

  // Participant-level priors
  w_0_pr ~ normal(0, 1);

  // Cutpoint priors
  cutpoints ~ normal(0, 1);

  // Ordered logistic likelihood with session as a fixed effect
  for (t in 1:T) {
    real pred_mean = w_0[subj_ix[t]] + beta_session * session[t];
    Y[t] ~ ordered_logistic(pred_mean, cutpoints);
  }
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T) {
    real pred_mean = w_0[subj_ix[t]] + beta_session * session[t];
    log_lik[t] = ordered_logistic_lpmf(Y[t] | pred_mean, cutpoints);
  }
}

