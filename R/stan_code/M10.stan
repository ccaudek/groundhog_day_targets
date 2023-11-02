data {
  int<lower=1> N; // total number of observations
  vector[N] Y; // response variable
  vector[N] trial; // centered trial predictor
  // data for group-level effects of ID 1 (random intercepts and slopes for subjects)
  int<lower=1> N_subjects; // number of subjects
  array[N] int<lower=1> subject; // subject index for each observation
}

parameters {
  real beta_pre; // slope before the change point
  real beta_post; // slope after the change point
  real<lower=0> sigma; // residual SD
  vector[N_subjects] subject_intercept; // random intercepts for subjects
  vector[N_subjects] subject_slope_pre; // random slopes before the change point
  vector[N_subjects] subject_slope_post; // random slopes after the change point
  real<lower=0> sigma_subject_intercept; // SD of subject intercepts
  real<lower=0> sigma_subject_slope_pre; // SD of subject slopes before the change point
  real<lower=0> sigma_subject_slope_post; // SD of subject slopes after the change point
  // Correlation between intercepts and slopes
  cholesky_factor_corr[3] L_subject;
}

model {
  // Priors
  beta_pre ~ normal(0, 10);
  beta_post ~ normal(0, 10);
  sigma ~ cauchy(0, 2.5);
  subject_intercept ~ normal(0, sigma_subject_intercept);
  subject_slope_pre ~ normal(0, sigma_subject_slope_pre);
  subject_slope_post ~ normal(0, sigma_subject_slope_post);
  L_subject ~ lkj_corr_cholesky(2);
  // SDs
  sigma_subject_intercept ~ cauchy(0, 2.5);
  sigma_subject_slope_pre ~ cauchy(0, 2.5);
  sigma_subject_slope_post ~ cauchy(0, 2.5);

  // Likelihood
  for (n in 1:N) {
    real mu = subject_intercept[subject[n]];
    // Add the random slopes depending on the phase of the trial
    mu += (trial[n] < 0) ? (beta_pre + subject_slope_pre[subject[n]]) * trial[n] 
                         : (beta_post + subject_slope_post[subject[n]]) * trial[n];
    Y[n] ~ normal(mu, sigma);
  }
}

generated quantities {
  real intercept; // intercept at the change point (trial == 0)
  vector[N] y_pred; // predicted values for each observation
  intercept = mean(subject_intercept); // assuming mean intercept across subjects

  // Predictions for each observation
  for (n in 1:N) {
    y_pred[n] = subject_intercept[subject[n]] + 
                ((trial[n] < 0) ? 
                (beta_pre + subject_slope_pre[subject[n]]) * trial[n] : 
                (beta_post + subject_slope_post[subject[n]]) * trial[n]);
  }

  // Compute the covariance matrix from the Cholesky factor
  cov_matrix[3] cov_subject = multiply_lower_tri_self_transpose(L_subject);
  // Extract the variances and correlations
  real var_subject_intercept = square(sigma_subject_intercept);
  real var_subject_slope_pre = square(sigma_subject_slope_pre);
  real var_subject_slope_post = square(sigma_subject_slope_post);
  real corr_intercept_slope_pre = L_subject[2, 1] * sigma_subject_intercept * sigma_subject_slope_pre;
  real corr_intercept_slope_post = L_subject[3, 1] * sigma_subject_intercept * sigma_subject_slope_post;
  real corr_slope_pre_slope_post = L_subject[3, 2] * sigma_subject_slope_pre * sigma_subject_slope_post;
}

