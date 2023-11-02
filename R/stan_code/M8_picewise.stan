data {
  int<lower=0> N;
  vector[N] happiness;
  matrix[N, 6] X;
  // ... (rest of your data block)
}

parameters {
  real alpha;
  vector[6] beta;
  real<lower=0> sigma;
  real<lower=2> nu;  // new parameter for degrees of freedom
}

transformed parameters {
  vector[N] mu;
  // Comment out z_subject

  mu = alpha + X * beta;
  
  // Comment out the loop for random effects
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ cauchy(0, 2);
  nu ~ gamma(2, 0.1);  // You can choose an appropriate prior for nu

  for (n in 1:N) {
    target += student_t_lpdf(happiness[n] | nu, mu[n], sigma);
  }
}