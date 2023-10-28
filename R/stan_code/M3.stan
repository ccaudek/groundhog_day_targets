data {
  int N;  // Total number of observations
  int J;  // Number of subjects
  array[J] int num_sessions;  // Number of sessions for each subject
  array[N] real happiness;  // Observed happiness for each observation
  array[N] real outcome_history;  // Outcome history for each observation
}

parameters {
  real w0;
  real w1;
  real<lower=0> sigma;  // Residual standard deviation
}

model {
  // Priors
  w0 ~ normal(0, 1);
  w1 ~ normal(0, 1);
  sigma ~ normal(0, 1);

  // Likelihood
  int n = 1;
  for (j in 1:J) {
    for (s in 1:num_sessions[j]) {
      for (t in 1:30) {  // Assuming each session has 30 trials
        real pred_happiness = w0 + w1 * outcome_history[n];
        happiness[n] ~ normal(pred_happiness, sigma);
        n += 1;
      }
    }
  }
}

generated quantities {
  array[N] real pred_happiness;  // For N observations
  
  int n = 1;
  for (j in 1:J) {
    for (s in 1:num_sessions[j]) {
      for (t in 1:30) {
        pred_happiness[n] = w0 + w1 * outcome_history[n];
        n += 1;
      }
    }
  }
}
