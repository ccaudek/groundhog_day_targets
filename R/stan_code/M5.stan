data {
  int<lower=0> N; // Number of observations
  int<lower=1> S; // Number of subjects
  int<lower=1> E; // Number of unique subject-session combinations
  vector[N] weighted_outcome; // Changed to vector
  vector[N] weighted_stimulus; // Changed to vector
  vector[N] weighted_RPE; // Changed to vector
  vector[N] zmoodpre;
  vector[N] zcontrol;
  vector[N] trial;
  vector[N] happiness;
  array[N] int<lower=1, upper=S> subject; // Subject index for each observation
  array[N] int<lower=1, upper=E> session; // Session index for each observation
}

parameters {
  // Group-level parameters
  real mu_alpha;
  vector[6] mu_beta;
  real<lower=0> sigma_alpha;
  vector<lower=0>[6] sigma_beta;

  // Subject-level parameters (non-centered)
  vector[S] alpha_subject_raw;
  matrix[6, S] beta_subject_raw;

  // Session-level parameters (non-centered)
  vector[E] alpha_session_raw;
  matrix[6, E] beta_session_raw;

  real<lower=0, upper=1> gamma;  // Exponential decay factor
}

transformed parameters {
  // Subject-level parameters (centered)
  vector[S] alpha_subject = mu_alpha + alpha_subject_raw * sigma_alpha;
  matrix[6, S] beta_subject;

  // Session-level parameters (centered)
  vector[E] alpha_session = mu_alpha + alpha_session_raw * sigma_alpha;
  matrix[6, E] beta_session;

  for (j in 1:6) {
    beta_subject[j] = mu_beta[j] + beta_subject_raw[j] * sigma_beta[j];
    beta_session[j] = mu_beta[j] + beta_session_raw[j] * sigma_beta[j];
  }
}

model {
  // Group-level priors
  mu_alpha ~ normal(0, 1);
  mu_beta ~ normal(0, 1);
  sigma_alpha ~ cauchy(0, 2.5);
  sigma_beta ~ cauchy(0, 2.5);

  // Non-centered subject-level priors
  alpha_subject_raw ~ std_normal();
  for (j in 1:6) {
    beta_subject_raw[j] ~ std_normal();
  }

  // Non-centered session-level priors
  alpha_session_raw ~ std_normal();
  for (j in 1:6) {
    beta_session_raw[j] ~ std_normal();
  }

  // Vectorized Likelihood
  vector[N] mu_happiness;
  for (n in 1:N) {
    int subj_idx = subject[n];
    int sess_idx = session[n];
    real subj_alpha = alpha_subject[subj_idx];
    real sess_alpha = alpha_session[sess_idx];
    
    mu_happiness[n] = subj_alpha + sess_alpha + 
      beta_subject[1, subj_idx] * weighted_outcome[n] +
      beta_subject[2, subj_idx] * weighted_stimulus[n] +
      beta_subject[3, subj_idx] * weighted_RPE[n] +
      beta_subject[4, subj_idx] * zmoodpre[n] +
      beta_subject[5, subj_idx] * zcontrol[n] +
      beta_subject[6, subj_idx] * trial[n];
  }
  happiness ~ normal(mu_happiness, 1);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for each observation for LOO
  vector[N] y_pred;   // Predicted values for happiness

  for (n in 1:N) {
    int subj_idx = subject[n];
    int sess_idx = session[n];
    real subj_alpha = alpha_subject[subj_idx];
    real sess_alpha = alpha_session[sess_idx];
    
    // Compute the mean of the Gaussian for each observation
    real mu_happiness = subj_alpha + sess_alpha + 
      beta_subject[1, subj_idx] * weighted_outcome[n] +
      beta_subject[2, subj_idx] * weighted_stimulus[n] +
      beta_subject[3, subj_idx] * weighted_RPE[n] +
      beta_subject[4, subj_idx] * zmoodpre[n] +
      beta_subject[5, subj_idx] * zcontrol[n] +
      beta_subject[6, subj_idx] * trial[n];

    // Compute the log-likelihood for each observation
    log_lik[n] = normal_lpdf(happiness[n] | mu_happiness, 1);

    // Generate predicted values for happiness
    y_pred[n] = normal_rng(mu_happiness, 1);
  }
}

