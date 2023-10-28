

# Prepare the data (this is a mock example; you should replace this with your real data)
N <- 100 # total number of observations
J <- 3  # total number of participants

stan_data <- list(
  N = N,
  J = J,
  participant = sample(1:J, N, replace = TRUE),
  outcome = rnorm(N),
  stimulus = rnorm(N),
  RPE = rnorm(N),
  zmoodpre = rnorm(N),
  zcontrol = rnorm(N),
  trial = rnorm(N),
  happiness = runif(N)
)


mod <- cmdstan_model("R/stan_code/M3.stan")

# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 1000,
  iter_warmup = 500,
  chains = 4
)


# Assuming 'mod' is your compiled Stan model and 'stan_data' is your data list
# fit <- mod$vb(
fit <- mod$variational(
  data = stan_data,
  iter = 50000,  # Number of iterations
  tol_rel_obj = 0.01,  # Convergence tolerance; default is 0.01
  output_samples = 1000  # Number of posterior samples to draw from the approximate posterior
)


# Summary of the fit
fit$summary()

# Check diagnostics
fit$diagnose()

