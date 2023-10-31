data {
  int<lower=0> N; // Number of trials
  vector[N] feedback; // Feedback for each trial
  vector[N] mood; // Mood for each trial
}

parameters {
  real alpha; // Intercept
  real beta; // Coefficient for feedback
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);

  // Likelihood
  mood ~ normal(alpha + beta * feedback, 1);
}
