data {
  int<lower=0> N;  // Number of participants
  int<lower=0> T;  // Number of trials
  array[T] int<lower=1, upper=N> subj_ix;  // Participant index for each datapoint
  array[T] int<lower=1, upper=5> Y;  // Dependent variable: ordinal response (1 to 5)
}

parameters {
  real w_0_mu_pr; // Group-level mean of intercepts
  real<lower=0> w_0_sigma_pr; // Group-level SD of intercepts
  
  vector[N] w_0_pr;  // Participant-level intercepts
  ordered[4] cutpoints;  // Cutpoints for ordinal categories
}

transformed parameters {
  vector[N] w_0;

  for (p_ix in 1:N){
    w_0[p_ix] = w_0_mu_pr + w_0_sigma_pr * w_0_pr[p_ix];
  }
}

model {
  // Group-level priors
  w_0_mu_pr ~ normal(0, 1);
  w_0_sigma_pr ~ exponential(0.1);

  // Participant-level priors
  w_0_pr ~ normal(0, 1);

  // Cutpoint priors
  cutpoints ~ normal(0, 1);

  // Ordered logistic likelihood
  for (t in 1:T) {
    Y[t] ~ ordered_logistic(w_0[subj_ix[t]], cutpoints);
  }
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T) {
    log_lik[t] = ordered_logistic_lpmf(Y[t] | w_0[subj_ix[t]], cutpoints);
  }
}
