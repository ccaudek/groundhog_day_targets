data {
  int N;  // Number of observations
  int J;  // Number of subjects
  int S;  // Maximum number of sessions across all subjects
  array[N] int subj;  // Subject index for each observation
  array[N] int session;  // Session index for each observation
  array[N] int num_trials;  // Number of trials for each session
  array[J] int num_sessions;  // Number of sessions for each subject
  array[N] real happiness;  // Observed happiness
  array[N] real outcome_history;  // Outcome history for each observation
  array[N] real stimulus_history;  // Stimulus history for each observation
  array[N] real RPE_history;  // RPE history for each observation
  array[N] real zmoodpre;  // Pre-mood for each observation
  array[N] real zcontrol;  // Control for each observation
}

parameters {
  real mu_w0; real<lower=0> sigma_w0;
  real mu_w1; real<lower=0> sigma_w1;
  real mu_w2; real<lower=0> sigma_w2;
  real mu_w3; real<lower=0> sigma_w3;
  real mu_w4; real<lower=0> sigma_w4;
  real mu_w5; real<lower=0> sigma_w5;
  real mu_w6; real<lower=0> sigma_w6;
  real mu_gamma; real<lower=0> sigma_gamma;
  real<lower=0> sigma;  // Residual standard deviation

  // Variable-length arrays for session-level parameters
  array[J] vector[S] w0;
  array[J] vector[S] w1;
  array[J] vector[S] w2;
  array[J] vector[S] w3;
  array[J] vector[S] w4;
  array[J] vector[S] w5;
  array[J] vector[S] w6;
  array[J] vector<lower=0, upper=1>[S] gamma;
}

model {
  // Priors for group-level parameters
  mu_w0 ~ normal(0, 1); sigma_w0 ~ exponential(1);
  mu_w1 ~ normal(0, 1); sigma_w1 ~ exponential(1);
  mu_w2 ~ normal(0, 1); sigma_w2 ~ exponential(1);
  mu_w3 ~ normal(0, 1); sigma_w3 ~ exponential(1);
  mu_w4 ~ normal(0, 1); sigma_w4 ~ exponential(1);
  mu_w5 ~ normal(0, 1); sigma_w5 ~ exponential(1);
  mu_w6 ~ normal(0, 1); sigma_w6 ~ exponential(1);
  mu_gamma ~ normal(0.5, 0.1); sigma_gamma ~ exponential(1);
  sigma ~ normal(0, 1);
  
  int n = 1;  // Initialize index for observations

  for (j in 1:J) {
    for (s in 1:num_sessions[j]) {
      for (t in 1:num_trials[n]) {

        // Your existing code to update predicted happiness with exponential decay
        real pred_happiness = 0;
        pred_happiness += pow(gamma[j][s], t-1) * (
          w0[j][s] +
          w1[j][s] * outcome_history[n] +
          w2[j][s] * stimulus_history[n] +
          w3[j][s] * RPE_history[n] +
          w4[j][s] * zmoodpre[n] +
          w5[j][s] * zcontrol[n] +
          w6[j][s] * num_trials[n]
        );

        happiness[n] ~ normal(pred_happiness, sigma);
        
        if (n < N) {
          n += 1;  // Increment index for observations
        }
      }
    }
  }
}


// You can also add the `generated quantities` block if you need it





// generated quantities {
//   array[N] real pred_happiness;  // Array to hold predicted happiness for each observation
// 
//   for (n in 1:N) {
//     pred_happiness[n] = 0;
//     for (t in 1:T) {
//       if (t <= num_sessions[subj[n]]) {
//         pred_happiness[n] += pow(gamma[subj[n]][t], t-1) * (
//           w0[subj[n]][t] +
//           w1[subj[n]][t] * outcome_history[n, t] +
//           w2[subj[n]][t] * stimulus_history[n, t] +
//           w3[subj[n]][t] * RPE_history[n, t] +
//           w4[subj[n]][t] * zmoodpre[n] +
//           w5[subj[n]][t] * zcontrol[n] +
//           w6[subj[n]][t] * trial[n]
//         );
//       }
//     }
//   }
// }


