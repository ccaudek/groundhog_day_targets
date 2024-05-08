data {
  int<lower=1> N; // Total number of observations
  int<lower=1> S; // Total number of subjects
  int<lower=1> T; // Number of trials per session (constant across subjects and sessions)
  array[S] int<lower=1> n_sessions; // Number of sessions per subject
  vector[N] x; // Trial number
  vector[N] y; // Response variable
  array[N] int<lower=1> subject; // Subject indicator
  array[N] int<lower=1> session; // Session indicator
}
parameters {
  real mu_alpha; // Population-level mean intercept
  real mu_beta1; // Population-level mean slope for the first segment
  real mu_beta2; // Population-level mean slope for the second segment
  real<lower=0> sigma_alpha; // Standard deviation of subject-specific intercepts
  real<lower=0> sigma_beta1; // Standard deviation of subject-specific slopes for the first segment
  real<lower=0> sigma_beta2; // Standard deviation of subject-specific slopes for the second segment
  vector[S] alpha_raw; // Standardized random intercepts for each subject
  vector[S] beta1_raw; // Standardized random slopes for the first segment for each subject
  vector[S] beta2_raw; // Standardized random slopes for the second segment for each subject
  real<lower=0> sigma; // Residual standard deviation
}
transformed parameters {
  vector[N] mu; // Expected value of y
  vector[S] alpha = mu_alpha + sigma_alpha * alpha_raw; // Subject-specific intercepts
  vector[S] beta1 = mu_beta1 + sigma_beta1 * beta1_raw; // Subject-specific slopes for the first segment
  vector[S] beta2 = mu_beta2 + sigma_beta2 * beta2_raw; // Subject-specific slopes for the second segment
  
  for (n in 1:N) {
    if (x[n] <= 16) {
      mu[n] = alpha[subject[n]] + beta1[subject[n]] * x[n];
    } else {
      mu[n] = alpha[subject[n]] + beta1[subject[n]] * 16 + beta2[subject[n]] * (x[n] - 16);
    }
  }
}
model {
  // Priors
  mu_alpha ~ normal(0, 5);
  mu_beta1 ~ normal(0, 5);
  mu_beta2 ~ normal(0, 5);
  sigma_alpha ~ exponential(1);
  sigma_beta1 ~ exponential(1);
  sigma_beta2 ~ exponential(1);
  alpha_raw ~ std_normal();
  beta1_raw ~ std_normal();
  beta2_raw ~ std_normal();
  sigma ~ exponential(1);
  
  // Likelihood
  y ~ normal(mu, sigma);
}
