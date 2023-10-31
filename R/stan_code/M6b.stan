data {
  int<lower=0> N;  // Number of observations
  int<lower=1> S;  // Number of subjects
  int<lower=1> E;  // Number of unique subject-session combinations
  int<lower=1> T;  // Number of trials
  vector[N] RPE; 
  vector[N] outcome; 
  vector[N] stimulus;
  vector[N] zmoodpre;
  vector[N] zcontrol;
  vector[N] trial;
  vector[N] happiness;
  array[N] int<lower=1, upper=S> subject;  // Subject index for each observation
  array[N] int<lower=1, upper=E> session;  // Session index for each observation
}

parameters {
  real mu_alpha;  // Group-level intercept
  vector[7] mu_beta;  // Group-level slopes (added one more for the interaction term)
  real<lower=0> sigma_alpha;  // Group-level standard deviation for intercept
  vector<lower=0>[7] sigma_beta;  // Group-level standard deviation for slopes
  real<lower=0, upper=1> gamma;  // Decay parameter

  vector[S] alpha_subject_raw;  // Subject-level intercepts (non-centered)
  matrix[7, S] beta_subject_raw;  // Subject-level slopes (non-centered, 7 instead of 6)

  vector[E] alpha_session_raw;  // Session-level intercepts (non-centered)
  matrix[7, E] beta_session_raw;  // Session-level slopes (non-centered, 7 instead of 6)
}

transformed parameters {
  vector[S] alpha_subject = mu_alpha + alpha_subject_raw * sigma_alpha;  // Centered
  matrix[7, S] beta_subject;  // Centered subject-level slopes

  vector[E] alpha_session = mu_alpha + alpha_session_raw * sigma_alpha;  // Centered
  matrix[7, E] beta_session;  // Centered session-level slopes

  for (j in 1:7) {
    beta_subject[j] = mu_beta[j] + beta_subject_raw[j] * sigma_beta[j];
    beta_session[j] = mu_beta[j] + beta_session_raw[j] * sigma_beta[j];
  }
}

model {
  mu_alpha ~ normal(0, 1);
  mu_beta ~ normal(0, 1);
  sigma_alpha ~ cauchy(0, 2.5);
  sigma_beta ~ cauchy(0, 2.5);
  gamma ~ beta(1, 1);  // Prior for gamma
  
  alpha_subject_raw ~ std_normal();
  for (j in 1:7) {
    beta_subject_raw[j] ~ std_normal();
  }

  alpha_session_raw ~ std_normal();
  for (j in 1:7) {
    beta_session_raw[j] ~ std_normal();
  }

  vector[N] mu_happiness;
  for (n in 1:N) {
    int subj_idx = subject[n];
    int sess_idx = session[n];
    real subj_alpha = alpha_subject[subj_idx];
    real sess_alpha = alpha_session[sess_idx];

    // Compute weighted values with exponential decay
    real weighted_outcome_local = pow(gamma, trial[n]-1) * outcome[n];
    real weighted_stimulus_local = pow(gamma, trial[n]-1) * stimulus[n];
    real weighted_RPE_local = pow(gamma, trial[n]-1) * RPE[n];

    // Compute the mean of the Gaussian for each observation
    mu_happiness[n] = subj_alpha + sess_alpha + 
      beta_subject[1, subj_idx] * weighted_outcome_local +
      beta_subject[2, subj_idx] * weighted_stimulus_local +
      beta_subject[3, subj_idx] * weighted_RPE_local +
      beta_subject[4, subj_idx] * zmoodpre[n] +
      beta_subject[5, subj_idx] * zcontrol[n] +
      beta_subject[6, subj_idx] * trial[n] +
      beta_subject[7, subj_idx] * trial[n] * weighted_RPE_local;  // Interaction term
  }
  happiness ~ normal(mu_happiness, 1);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for LOO
  vector[N] y_pred;   // Predicted values

  for (n in 1:N) {
    int subj_idx = subject[n];
    int sess_idx = session[n];
    real subj_alpha = alpha_subject[subj_idx];
    real sess_alpha = alpha_session[sess_idx];

    // Compute weighted values with exponential decay
    real weighted_outcome_local = pow(gamma, trial[n]-1) * outcome[n];
    real weighted_stimulus_local = pow(gamma, trial[n]-1) * stimulus[n];
    real weighted_RPE_local = pow(gamma, trial[n]-1) * RPE[n];

    // Compute the mean of the Gaussian for each observation
    real mu_happiness = subj_alpha + sess_alpha + 
      beta_subject[1, subj_idx] * weighted_outcome_local +
      beta_subject[2, subj_idx] * weighted_stimulus_local +
      beta_subject[3, subj_idx] * weighted_RPE_local +
      beta_subject[4, subj_idx] * zmoodpre[n] +
      beta_subject[5, subj_idx] * zcontrol[n] +
      beta_subject[6, subj_idx] * trial[n] +
      beta_subject[7, subj_idx] * trial[n] * weighted_RPE_local;  // Interaction term

    log_lik[n] = normal_lpdf(happiness[n] | mu_happiness, 1);
    y_pred[n] = normal_rng(mu_happiness, 1);
  }
}
