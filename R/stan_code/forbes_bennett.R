
library(cmdstanr)
library(posterior)
library(bayesplot)
library(rstan)
library(loo)
color_scheme_set("brightblue")

cmdstan_version()


# M0 multiple subjects, only one session ---------------------------------------

# Preparing the data for Stan in a more robust and efficient manner
prepare_stan_data <- function(df) {
  
  # Filter and select columns, sort by 'user_id' and 'trial'
  prepared_df <- df %>%
    filter(ema_number == 1) %>%
    select(user_id, trial, instant_mood) %>%
    arrange(user_id, trial)
  
  # Convert 'user_id' to a numerical index for Stan
  prepared_df$subj_ix <- as.integer(as.factor(prepared_df$user_id))
  
  # Make trial numbers zero-based
  prepared_df$trial_num <- prepared_df$trial - 1
  
  # Number of unique subjects and trials per subject
  N <- n_distinct(prepared_df$user_id)
  T_per_subject <- max(prepared_df$trial_num) + 1
  T <- N * T_per_subject
  
  # Create the list to be passed to cmdstan
  stan_data <- list(
    N = N,
    T = T,
    subj_ix = prepared_df$subj_ix,
    trial_num = prepared_df$trial_num,
    Y = prepared_df$instant_mood
  )
  
  return(stan_data)
}


# Identify unique user_ids
unique_user_ids <- unique(prl_df$user_id)

# Randomly sample 10 user_ids
random_user_ids <- sample(unique_user_ids, 10)

# Filter the original data frame to include only these 10 user_ids
filtered_df <- prl_df %>% 
  filter(user_id %in% random_user_ids)

stan_data <- prepare_stan_data(filtered_df)

mod <- cmdstan_model("R/stan_code/M0.stan")

# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  parallel_chains = 4,
  chains = 4,
  refresh = 500
)

fit$diagnostic_summary()

# In order to use this method you must compute and save the pointwise 
# log-likelihood in your Stan program.
loo_m0 <- fit$loo(cores = 4)
print(loo_m0)

fit$summary()

fit_vb <- mod$variational(
  data = stan_data,
  seed = 123
)

fit$diagnostic_summary()
fit$summary()

summary_fit <- fit$summary()
# Display the column names of the summary data frame
colnames(summary_fit)
summary_fit$variable


# M1_fixed ---------------------------------------------------------------------
#' Multiple subjects, multiple sessions. Session as fixed-effect.

prepare_stan_data_with_session <- function(df) {
  
  # Filter and select columns, sort by 'user_id' and 'trial'
  prepared_df <- df %>%
    select(user_id, trial, instant_mood, ema_number) %>%  
    arrange(user_id, ema_number, trial)  
  
  # Convert 'user_id' to a numerical index for Stan
  prepared_df$subj_ix <- as.integer(as.factor(prepared_df$user_id))
  
  # Make trial numbers zero-based
  prepared_df$trial_num <- prepared_df$trial - 1
  
  # Number of unique subjects and trials per subject
  N <- n_distinct(prepared_df$user_id)
  T_per_subject <- max(prepared_df$trial_num) + 1
  T <- nrow(prepared_df)
  
  # Create a data frame that contains session-level information
  session_df <- prepared_df %>%
    select(user_id, subj_ix, ema_number) %>%
    distinct() %>%
    arrange(user_id, ema_number)
  
  # Total number of sessions across all subjects
  S <- nrow(session_df)
  
  # Subject index for each session
  session_subj <- session_df$subj_ix
  
  # Create the list to be passed to cmdstan
  stan_data <- list(
    N = N,
    T = T,
    subj_ix = prepared_df$subj_ix,
    trial_num = prepared_df$trial_num,
    Y = prepared_df$instant_mood,
    session = prepared_df$ema_number,
    S = S,  # Add total number of sessions
    session_subj = session_subj  # Add subject index for each session
  )
  
  return(stan_data)
}


set.seed(42)
# Identify unique user_ids
unique_user_ids <- unique(prl_df$user_id)
# Randomly sample 10 user_ids
random_user_ids <- sample(unique_user_ids, 3)
stan_data <- prepare_stan_data_with_session(filtered_df)

mod1f <- cmdstan_model("R/stan_code/M1_fixed.stan")

fit_m1_fixed_vb <- mod1f$variational(
  data = stan_data,
  seed = 12345
)

fit1f_vb$diagnostic_summary()
fit1f_vb$summary()


# Fit the model to the data
fit1f <- mod1f$sample(
  data = stan_data,
  iter_sampling = 200,
  iter_warmup = 100,
  parallel_chains = 4,
  chains = 4
)

fit1f$diagnostic_summary()
fit1f$summary()

mcmc_hist(fit1f$draws("beta_session"))
mcmc_hist(fit1f$draws("beta_session"))

loo_m1f <- fit1f$loo(cores = 4)
print(loo_m1f)


# M1_random --------------------------------------------------------------------

mod1r <- cmdstan_model("R/stan_code/M1_random.stan")

fit_m1_random_vb <- mod1r$variational(
  data = stan_data,
  seed = 12345
)

fit_m1_random_vb$summary()


# Fit the model to the data
fit1r <- mod1r$sample(
  data = stan_data,
  iter_sampling = 200,
  iter_warmup = 100,
  parallel_chains = 4,
  chains = 4
)

fit1r$diagnostic_summary()
fit1r$summary()

mcmc_hist(fit1r$draws("beta_session"))

loo_m1r <- fit1r$loo(cores = 4)
print(loo_m1r)

loo_compare(loo_m1f, loo_m1r)


# M2r_rpe ----------------------------------------------------------------------

#' Get the necessary information from the prl_df and the 
#' params_happiness_df data frames.
from_prl_df <- prl_df |> 
  dplyr::select(
    user_id, is_reversal, ema_number, trial,
    instant_mood, feedback, is_target_chosen
  )

params_happiness_df <- params_happiness_df |> 
  mutate(
    is_reversal = ifelse(is_reversal == 1, "yes", "no")
  )

all_data_df <- left_join(
  params_happiness_df, from_prl_df, 
  by = c("user_id", "is_reversal", "ema_number", "trial")
)


prepare_stan_data <- function(df) {
  
  # Filter and select columns, sort by 'user_id' and 'trial'
  prepared_df <- df %>%
    select(user_id, ema_number, trial, instant_mood, RPE) %>%  
    arrange(user_id, ema_number, trial)  
  
  # Convert 'user_id' to a numerical index for Stan
  prepared_df$subj_ix <- as.integer(as.factor(prepared_df$user_id))
  
  # Number of unique subjects and trials per subject
  N <- n_distinct(prepared_df$user_id)
  T_per_subject <- max(prepared_df$trial)
  T <- nrow(prepared_df)
  
  # Create a data frame that contains session-level information
  session_df <- prepared_df %>%
    select(user_id, subj_ix, ema_number) %>%
    distinct() %>%
    arrange(user_id, ema_number)
  
  # Total number of sessions across all subjects
  S <- nrow(session_df)
  
  # Subject index for each session
  session_subj <- session_df$subj_ix
  
  # Create the list to be passed to cmdstan
  stan_data <- list(
    N = N,
    T = T,
    subj_ix = prepared_df$subj_ix,
    Y = prepared_df$instant_mood,
    session = prepared_df$ema_number,
    S = S,
    session_subj = session_subj,
    RPE = prepared_df$RPE
  )
  
  return(stan_data)
}


# Identify unique user_ids
unique_user_ids <- unique(all_data_df$user_id)
# Randomly sample 10 user_ids
random_user_ids <- sample(unique_user_ids, 3)
# Filter the original data frame to include only these 10 user_ids
filtered_df <- all_data_df %>% 
  filter(user_id %in% random_user_ids)

stan_data <- prepare_stan_data(filtered_df)

mod2r <- cmdstan_model("R/stan_code/M2r_rpe.stan")

fit_mod2r_vb <- mod2r$variational(
  data = stan_data,
  seed = 42
)

fit_m1_random_vb$summary()



# Fit the model to the data
fit2r <- mod2r$sample(
  data = stan_data,
  iter_sampling = 500,
  iter_warmup = 200,
  parallel_chains = 4,
  chains = 4
)

fit2r$diagnostic_summary()
fit2r$summary()

summary_out <- fit2r$summary()
summary_out$variable

mcmc_hist(fit2r$draws("beta_RPE"))



















