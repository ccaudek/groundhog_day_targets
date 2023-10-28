

library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

cmdstan_version()




tar_load(params_happiness_df)

user_id_list <- unique(params_happiness_df$user_id)



subset_df <- 
  params_happiness_df[params_happiness_df$user_id %in% user_id_list[1:3], ]

# Sort the data frame
sorted_df <- subset_df %>%
  arrange(user_id, ema_number, trial)


# Number of observations
N <- nrow(sorted_df)

# Number of subjects
J <- length(unique(sorted_df$user_id))

# Maximum number of sessions
S <- max(sorted_df$ema_number)

# Subject index for each observation
subj <- as.integer(factor(sorted_df$user_id))

# Session index for each observation
session <- sorted_df$ema_number

# Extract other variables
happiness <- sorted_df$happiness
outcome_history <- as.matrix(sorted_df %>% select(outcome))
stimulus_history <- as.matrix(sorted_df %>% select(stimulus))
RPE_history <- as.matrix(sorted_df %>% select(RPE))
zmoodpre <- sorted_df$zmoodpre
zcontrol <- sorted_df$zcontrol
trial <- sorted_df$trial

# Number of sessions for each subject
num_sessions_table <- table(sorted_df$user_id)
num_sessions <- as.integer(num_sessions_table) / 30  # Assuming 30 trials per session
names(num_sessions) <- names(num_sessions_table)

# Number of trials for each subject (this would require your specific logic)
# For example, if each session has a fixed number of trials (let's say 30),
# then num_trials for each subject would be 30 times the number of sessions for that subject.
num_trials <- num_sessions * 30  # Adjust this calculation as needed

# Ensure that the order of 'num_sessions' and 'num_trials' matches the order of subjects in 'subj'
ordered_num_sessions <- num_sessions[as.character(unique(sorted_df$user_id))]
ordered_num_trials <- num_trials[as.character(unique(sorted_df$user_id))]

# Convert one-column matrices to vectors
outcome_history_vector <- as.vector(sorted_df$outcome)
stimulus_history_vector <- as.vector(sorted_df$stimulus)
RPE_history_vector <- as.vector(sorted_df$RPE)

# Create the data list with the added 'num_sessions', 'num_trials', and 'S'
stan_data <- list(
  N = N,
  J = J,
  S = S,
  subj = subj,
  session = session,
  num_sessions = ordered_num_sessions,
  num_trials = ordered_num_trials,
  happiness = happiness,
  outcome_history = outcome_history_vector,
  stimulus_history = stimulus_history_vector,
  RPE_history = RPE_history_vector,
  zmoodpre = zmoodpre,
  zcontrol = zcontrol,
  trial = trial
)

stan_data$num_trials <- rep(30, 480)  # Replace 30 with the correct value

# Validate that the length is now 480
print(length(stan_data$num_trials))


# Now you can use stan_data as the data argument when calling the Stan function.


mod <- cmdstan_model("R/stan_code/M5.stan")


# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 100,
  iter_warmup = 50,
  chains = 4
)

# Summary of the fit
fit$summary()


fit_mle <- mod$optimize(data = stan_data, seed = 42)

fit_mle$print() # includes lp__ (log prob calculated by Stan program)
fit_mle$print("mu_w0")
fit_mle$print("mu_w1")
fit_mle$print("mu_w2")
fit_mle$print("mu_w3")
fit_mle$print("mu_w4")
fit_mle$print("mu_w5")
fit_mle$print("mu_w6")
fit_mle$print("mu_gamma")



fit_vb <- mod$variational(
  data = stan_data,
  seed = 42,
  iter = 40000
)

fit_vb$print("mu_w0")
fit_vb$print("mu_w1")
fit_vb$print("mu_w2")
fit_vb$print("mu_w3")
fit_vb$print("mu_w4")
fit_vb$print("mu_w5")
fit_vb$print("mu_w6")
fit_vb$print("mu_gamma")




##### This works! ###############

subset_df <- 
  params_happiness_df[params_happiness_df$user_id %in% user_id_list[1:3], ]

# Sort the data frame
sorted_df <- subset_df %>%
  arrange(user_id, ema_number, trial)

unique_users <- unique(sorted_df$user_id)
J <- length(unique_users)

# Calculate the number of sessions for each subject
num_sessions <- sapply(unique_users, function(u) {
  max(sorted_df[sorted_df$user_id == u, 'ema_number'])
})

# Prepare the stan_data list
stan_data <- list(
  N = nrow(sorted_df),
  J = J,
  num_sessions = num_sessions,
  happiness = sorted_df$happiness,
  outcome_history = sorted_df$outcome
)

mod <- cmdstan_model("R/stan_code/M3.stan")


# Fit the model to the data
fit_simplified <- mod$sample(
  data = stan_data,
  iter_sampling = 1000,
  iter_warmup = 500,
  chains = 4
)

# Check the dimensions of the posterior samples
posterior_samples_simplified <- fit_simplified$draws()
dim(posterior_samples_simplified)

fit_simplified$summary()


# Get the summary tibble
summary_tibble <- fit_simplified$summary()

# Filter for rows where the variable starts with 'pred_happiness'
pred_happiness_rows <- summary_tibble %>% 
  filter(grepl("^pred_happiness", variable))

# Extract the mean of the posterior distributions for these rows
pred_happiness_means <- pred_happiness_rows$mean
length(pred_happiness_means)

cor(
  stan_data$happiness,
  pred_happiness_means
)

# Get posterior draws
draws <- fit$draws()
print(draws)
as_draws_df(draws)
mcmc_hist(fit$draws("w1"))

fit$cmdstan_diagnose()

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

mcmc_intervals(stanfit, pars = c("w0", "w1"))
mcmc_areas(
  stanfit, 
  pars = c("w0", "w1"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "mean"
)

# Step 1: Extract Posterior Samples
posterior_samples <- fit$draws() # Assuming fit is your cmdstanr fit object
w0_samples <- posterior_samples[, , "w0"]
w1_samples <- posterior_samples[, , "w1"]
sigma_samples <- posterior_samples[, , "sigma"]

# Step 2: Generate yrep
y <- stan_data$happiness
x <- stan_data$outcome_history
N <- length(y) # Assuming y is your observed outcome data
n_samples <- dim(posterior_samples)[1] * dim(posterior_samples)[2] # Total number of posterior samples

# Initialize yrep matrix
yrep <- matrix(0, n_samples, N)

for (i in 1:n_samples) {
  w0 <- w0_samples[i]
  w1 <- w1_samples[i]
  sigma <- sigma_samples[i]
  
  # Generate posterior predictive sample for this set of parameters
  yrep[i, ] <- rnorm(N, w0 + w1 * x, sigma) # Replace 'x' with your predictor variable
}

# Step 3: Run pp_check
pp_check(y, yrep[1:50, ], ppc_dens_overlay)




