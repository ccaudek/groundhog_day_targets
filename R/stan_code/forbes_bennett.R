
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
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


# prepare_stan_data <- function(df) {
#   
#   # Filter and select columns, sort by 'user_id' and 'trial'
#   prepared_df <- df %>%
#     select(user_id, ema_number, trial, instant_mood, RPE) %>%  
#     arrange(user_id, ema_number, trial)  
#   
#   # Convert 'user_id' to a numerical index for Stan
#   prepared_df$subj_ix <- as.integer(as.factor(prepared_df$user_id))
#   
#   # Number of unique subjects and trials per subject
#   N <- n_distinct(prepared_df$user_id)
#   T_per_subject <- max(prepared_df$trial)
#   T <- nrow(prepared_df)
#   
#   # Create a data frame that contains session-level information
#   session_df <- prepared_df %>%
#     select(user_id, subj_ix, ema_number) %>%
#     distinct() %>%
#     arrange(user_id, ema_number)
#   
#   # Total number of sessions across all subjects
#   S <- nrow(session_df)
#   
#   # Subject index for each session
#   session_subj <- session_df$subj_ix
#   
#   # Create the list to be passed to cmdstan
#   stan_data <- list(
#     N = N,
#     T = T,
#     subj_ix = prepared_df$subj_ix,
#     Y = prepared_df$instant_mood,
#     session = prepared_df$ema_number,
#     S = S,
#     session_subj = session_subj,
#     RPE = prepared_df$RPE
#   )
#   
#   return(stan_data)
# }

prepare_stan_data <- function(df) {

  # Filter and select columns, sort by 'user_id' and 'trial'
  prepared_df <- df %>%
    select(user_id, ema_number, trial, instant_mood, RPE, outcome, stimulus, zmoodpre, zcontrol) %>%
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
  
  # Create a numerical index for ema_number
  prepared_df$session_ix <- as.integer(as.factor(prepared_df$ema_number))
  

  # Create the list to be passed to cmdstan
  stan_data <- list(
    N = N,
    T = T,
    subj_ix = prepared_df$subj_ix,
    Y = prepared_df$instant_mood,
    session = prepared_df$ema_number,
    S = S,
    session_subj = session_subj,
    RPE = prepared_df$RPE,
    outcome = prepared_df$outcome,
    stimulus = prepared_df$stimulus,
    zmoodpre = prepared_df$zmoodpre,
    zcontrol = prepared_df$zcontrol,
    trial = prepared_df$trial
  )

  return(stan_data)
} #### good for M3

##### This is for M4
# prepare_stan_data <- function(df) {
#   
#   # Filter and select columns, sort by 'user_id', 'ema_number', and 'trial'
#   prepared_df <- df %>%
#     select(user_id, ema_number, trial, happiness, RPE, outcome, stimulus, 
#            zmoodpre, zcontrol) %>%
#     arrange(user_id, ema_number, trial)
#   
#   # Convert 'user_id' to a numerical index for Stan
#   prepared_df$subj_ix <- as.integer(as.factor(prepared_df$user_id))
#   
#   # Number of unique subjects and total trials
#   N <- n_distinct(prepared_df$user_id)
#   T <- nrow(prepared_df)
#   
#   # Create a data frame that contains session-level information
#   session_df <- prepared_df %>%
#     select(user_id, subj_ix, ema_number) %>%
#     distinct() %>%
#     arrange(user_id, ema_number)
#   
#   # Total number of sessions across all subjects
#   S <- nrow(session_df)
#   
#   # Subject index for each session
#   session_subj <- session_df$subj_ix
#   
#   # Create the list to be passed to CmdStan
#   stan_data <- list(
#     N = N,
#     T = T,
#     subj_ix = prepared_df$subj_ix,
#     Y = prepared_df$happiness, 
#     session = prepared_df$ema_number,
#     S = S,
#     session_subj = session_subj,
#     RPE = prepared_df$RPE,
#     outcome = prepared_df$outcome,
#     stimulus = prepared_df$stimulus,
#     zmoodpre = prepared_df$zmoodpre,
#     zcontrol = prepared_df$zcontrol,
#     trial = prepared_df$trial
#   )
#   
#   return(stan_data)
# }




# Identify unique user_ids
unique_user_ids <- unique(all_data_df$user_id)
# Randomly sample 10 user_ids
random_user_ids <- sample(unique_user_ids, 2)
# Filter the original data frame to include only these 10 user_ids
filtered_df <- all_data_df %>% 
  filter(user_id %in% random_user_ids)

data_for_stan <- prepare_stan_data(filtered_df)
# stan_data <- prepare_stan_data(all_data_df)


# mod2r$format(
#   canonicalize = list("deprecations"),
#   overwrite_file = TRUE,
#   backup = FALSE
# )
# mod2r$print()

mod2r <- cmdstan_model(
  "R/stan_code/M2r_rpe_v2.stan", 
  stanc_options = list("O1"),
  force_recompile = TRUE
)

fit_mod2r_vb <- mod2r$variational(
  data = data_for_stan,
  seed = 42
)

fit_mod2r_vb$summary()

log_lik_vb <- fit_mod2r_vb$draws("log_lik")  

# Compute LOO
loo_fit_mod2r_vb <- loo(log_lik_vb)
print(loo_fit_mod2r_vb)

mcmc_hist(fit_mod2r_vb$draws("beta_RPE"))




















# oggi

sample_mod <- cmdstan_model(
  "R/stan_code/M3.stan", 
  stanc_options = list("O1"),
  force_recompile = TRUE
)


# Identify unique user_ids
set.seed(123)
unique_user_ids <- unique(all_data_df$user_id)
# Randomly sample 2 user_ids
random_user_ids <- sample(unique_user_ids, 50)
random_user_ids
# Filter the original data frame to include only these 10 user_ids
filtered_df <- all_data_df %>% 
  filter(user_id %in% random_user_ids)

data_for_stan <- prepare_stan_data(filtered_df)


# Check with VI whether the model converges
sampled_vb <- sample_mod$variational(
  data = data_for_stan,
  seed = 42
)

sampled_vb$summary()

log_lik_vb <- sampled_vb$draws("log_lik")  

# Compute LOO
loo_vb <- loo(log_lik_vb)
print(loo_vb)

summary_mod <- sampled_vb$summary()
summary_mod$variable
summary(summary_mod$rhat)




mcmc_hist(sampled_vb$draws("beta_RPE"))
mcmc_hist(sampled_vb$draws("beta_outcome"))
mcmc_hist(sampled_vb$draws("beta_stimulus"))
mcmc_hist(sampled_vb$draws("beta_zmoodpre"))
mcmc_hist(sampled_vb$draws("beta_zcontrol"))
mcmc_hist(sampled_vb$draws("beta_trial"))

# Sample with MCMC
sampled <- sample_mod$sample(
  data = data_for_stan,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 300,
  iter_warmup = 100,
  refresh = 20
)

sampled$diagnostic_summary()
sampled$summary()

# Compute LOO
loo_m3 <- sampled$loo(cores = 4)
print(loo_m3)

summary_mod <- sampled$summary()
summary_mod$variable
summary(summary_mod$rhat)

mcmc_hist(sampled$draws("beta_RPE"))
mcmc_trace(sampled$draws("beta_RPE"))
mcmc_hist(sampled$draws("beta_outcome"))
mcmc_hist(sampled$draws("beta_stimulus"))
mcmc_hist(sampled$draws("beta_zmoodpre"))
mcmc_hist(sampled$draws("beta_zcontrol"))
mcmc_hist(sampled$draws("beta_trial"))


# Extract the Samples
Y_pred_draws <- sampled_vb$draws(variables = "Y_pred")
# Convert to an array
Y_pred_array <- as.array(Y_pred_draws)

# Compute the mode for each trial
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) mlv(as.vector(x)))
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) median(as.vector(x)))
Y_pred_mode <- apply(Y_pred_array, 3, function(x) mean(as.vector(x), trim=0.1))

Y_pred_mode <- colMeans(Y_pred_array)

# Assume Y_observed is a vector of observed values
comparison <- data.frame(
  Y_observed = data_for_stan$Y, 
  Y_pred = Y_pred_mode
)

cor(
  comparison$Y_observed, comparison$Y_pred
  )


filtered_df$happiness_hat <- Y_pred_mode


# Calculate mean and standard error of the mean (SEM)
for_plot_df <- filtered_df |> 
  group_by(is_reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(happiness_hat, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(happiness_hat, na.rm = TRUE) / sqrt(n())
  ) |> 
  ungroup()
for_plot_df$is_reversal <- as.factor(for_plot_df$is_reversal)

# Create the plot
for_plot_df |> 
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = viridis(3)[1], linewidth = 2) +  # Line for h
  geom_ribbon(aes(ymin = h - h_sem, ymax = h + h_sem), alpha = 0.2) +  # Ribbon for h
  geom_line(aes(y=h_hat), color = viridis(3)[2], linewidth = 1.0) +  # Line for h_hat
  #geom_ribbon(aes(ymin = h_hat - h_hat_sem, ymax = h_hat + h_hat_sem), alpha = 0.2, fill = "red") +  # Ribbon for h_hat
  facet_wrap(~ is_reversal) +
  theme_default()  # Apply the bayesplot theme






# fit_map <- mod2r$optimize(data = stan_data, seed = 42)

# Load necessary libraries
library(dplyr)
library(tidyr)

# Assuming fit_map$summary() is stored in a variable called summary_data
summary_data <- fit_map$summary()

# Filter out the rows corresponding to the individual-level parameters
individual_params <- summary_data %>%
  filter(grepl("^w_0_pr\\[", variable)) %>%
  separate(variable, into = c("parameter", "subject"), sep = "\\[|\\]", remove = FALSE) %>%
  mutate(subject = as.integer(subject)) %>%
  select(subject, estimate)

# Reshape the data from long to wide format
individual_params_wide <- individual_params %>%
  group_by(subject) %>%
  mutate(row_num = row_number()) %>%
  pivot_wider(names_from = row_num, values_from = estimate, names_prefix = "w_0_pr_")

# Now individual_params_wide is a data frame with one row per subject and one column per parameter
print(individual_params_wide)








# fit_laplace <- mod2r$laplace(
#   mode = fit_map, 
#   draws = 100, 
#   data = stan_data, 
#   seed = 123, 
#   refresh = 10
# )


# DONE UNTIL HERE!

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
mcmc_hist(fit2r$draws("gamma"))

mcmc_hist(fit2r$draws("b_session[1]"))


fit2r_vb <- mod2r$variational(
  data = stan_data, 
  output_samples = 5000, 
  algorithm = "fullrank"
  )

fit2r_vb$summary()
mcmc_hist(fit2r_vb$draws("beta_RPE"))




#############################


##########

new_data <- params_happiness_df %>% 
  unnest_wider(mle_params, names_sep = "_")

library(dplyr)
library(tidyr)

# 1. Filter Columns
filtered_df <- params_happiness_df %>%
  select(is_reversal, happiness, mle_params, ema_number, user_id, best_alpha, mood_pre, mood_post, control)

# 2. Unnest mle_params into separate columns
unnested_df <- filtered_df %>%
  unnest_wider(mle_params, names_sep = "_")

# 3. Group by user_id and ema_number, then summarize
final_df <- unnested_df %>%
  group_by(user_id, ema_number) %>%
  summarise(across(
    where(is.numeric),
    ~mean(.x, na.rm = TRUE)  # Updated syntax
  ),
  is_reversal = first(is_reversal),
  .groups = 'drop'  # This drops the grouping, making it a regular data frame
  )

# Show the final data frame
glimpse(final_df)

library(lme4)

# Run the linear mixed-effects model
model <- lmer(happiness ~ mle_params_1 + mle_params_2 + mle_params_3 + mle_params_4 +
                mle_params_5 + mle_params_6 + mle_params_7 + mle_params_8 + 
                is_reversal + (1|user_id) + (1|ema_number), 
              data = final_df)

# Summary of the model to check the importance of each predictor
summary(model)

m <- brm(
  happiness ~ mle_params_1 + mle_params_2 + mle_params_3 + mle_params_4 +
    mle_params_5 + mle_params_6 + mle_params_7 + mle_params_8 + 
    is_reversal + (1|user_id) + (1|ema_number), 
  backend = "cmdstanr",
  family = student(),
  data = final_df
)

pp_check(m) + xlim(-3, 3)
bayes_R2(m)
summary(m)




# Self-compassion

tar_load(params_happiness_df)

new_data <- params_happiness_df %>% 
  unnest_wider(mle_params, names_sep = "_")


scs_scores <- rio::import(
  here::here("data", "prep", "quest_scales", "scs_scores.csv")
)
scs_scores$user_id <- as.numeric(scs_scores$user_id)

temp <- left_join(
  new_data, scs_scores, by = "user_id"
)

# Add happiness' model parameters names
names_cols <- colnames(temp)

weights_names <- c(
  "w0", "w_outcome", "w_stim", "w_rpe", "w_moodpre", "w_control",
  "w_ntrial", "w_gamma"
)

names_cols[12:19] <- weights_names
colnames(temp) <- names_cols

no_reversal_df <- temp |> 
  dplyr::filter(reversal == 0)

no_reversal_df$scs <- as.vector(scale(no_reversal_df$scs_total_score))

mod <- brm(
  happiness ~ 1 + scs * (w0 + w_outcome + w_stim + w_rpe + 
    w_moodpre + w_control + w_ntrial) + 
    (1 + w0 + w_outcome + w_stim + w_rpe + 
       w_moodpre + w_control + w_ntrial | user_id) + (1 | ema_number),
  family = student(),
  data = no_reversal_df,
  algorithm = "meanfield",
  iter = 40000
)
pp_check(mod) + xlim(-3, 3)

summary(mod)
conditional_effects(mod, "w0:scs")
conditional_effects(mod, "w_outcome:scs")
conditional_effects(mod, "w_stim:scs")
conditional_effects(mod, "w_moodpre:scs")
conditional_effects(mod, "w_control:scs")
conditional_effects(mod, "w_ntrial:scs")

bayes_R2(mod)


reversal_df <- temp |> 
  dplyr::filter(reversal == 1)

reversal_df$scs <- as.vector(scale(reversal_df$scs_total_score))


mod2 <- brm(
  happiness ~ 1 + scs * (w0 + w_outcome + w_stim + w_rpe + 
                           w_moodpre + w_control + w_ntrial) + 
    (1 + w0 + w_outcome + w_stim + w_rpe + 
       w_moodpre + w_control + w_ntrial | user_id) + (1 | ema_number),
  family = student(),
  data = reversal_df,
  algorithm = "meanfield",
  iter = 40000
)
pp_check(mod2) + xlim(-3, 3)

summary(mod2)
conditional_effects(mod2, "w_outcome:scs")
conditional_effects(mod2, "w_stim:scs")
conditional_effects(mod2, "w_rpe:scs")
conditional_effects(mod2, "w_moodpre:scs")
conditional_effects(mod2, "w_control:scs")
conditional_effects(mod2, "w_ntrial:scs")

bayes_R2(mod2)

temp$scs <- as.vector(scale(temp$scs_total_score))

mod3 <- brm(
  happiness ~ 1 + reversal + scs * 
    (w0 + w_outcome + w_stim + w_rpe + w_moodpre + w_control + w_ntrial) + 
    (w0 + w_outcome + w_stim + w_rpe + w_moodpre + w_control + 
       w_ntrial | user_id) + (1 | ema_number),
  family = asym_laplace(),
  data = temp,
  backend = "cmdstanr"
  # algorithm = "meanfield",
  #iter = 40000
)
pp_check(mod3) + xlim(-3, 3)
bayes_R2(mod3)

summary(mod3)
conditional_effects(mod3, "w_outcome:scs")
conditional_effects(mod3, "w_stim:scs")
conditional_effects(mod3, "w_rpe:scs")
conditional_effects(mod3, "w_moodpre:scs")
conditional_effects(mod3, "w_control:scs")
conditional_effects(mod3, "w_ntrial:scs")



################################################################################



prepare_stan_data <- function(df) {
  
  # Filter and select columns, sort by 'user_id' and 'trial'
  prepared_df <- df %>%
    select(user_id, ema_number, trial, happiness, RPE, outcome, stimulus, 
           zmoodpre, zcontrol) %>%
    arrange(user_id, ema_number, trial)
  
  # Convert 'user_id' to a numerical index for Stan
  prepared_df$subj_idx <- as.integer(as.factor(prepared_df$user_id))
  
  # Number of unique subjects and unique subject-session combinations
  N <- nrow(prepared_df)
  S <- length(unique(prepared_df$user_id))
  E <- nrow(unique(prepared_df[, c("user_id", "ema_number")]))
  T <- max(prepared_df$trial)
  
  # Create the list to be passed to cmdstan
  stan_data <- list(
    N = N,
    S = S,
    E = E,
    T = T,
    subject = prepared_df$subj_idx,
    session = prepared_df$ema_number,
    outcome = as.vector(prepared_df$outcome),  # Not weighted
    stimulus = as.vector(prepared_df$stimulus),  # Not weighted
    RPE = as.vector(prepared_df$RPE),  # Not weighted
    zmoodpre = prepared_df$zmoodpre,
    zcontrol = prepared_df$zcontrol,
    trial = prepared_df$trial,
    happiness = prepared_df$happiness
  )
  
  return(stan_data)
}



sample_mod <- cmdstan_model(
  "R/stan_code/M6b.stan", 
  stanc_options = list("O1"),
  force_recompile = TRUE
)

# Identify unique user_ids
set.seed(1)
unique_user_ids <- unique(all_data_df$user_id)
# Randomly sample 2 user_ids
random_user_ids <- sample(unique_user_ids, 3)
random_user_ids
# Filter the original data frame to include only these 10 user_ids
filtered_df <- all_data_df %>% 
  filter(user_id %in% random_user_ids)

# data_for_stan <- prepare_stan_data(filtered_df)
data_for_stan <- prepare_stan_data(all_data_df)


# Check with VI whether the model converges
sampled_vb <- sample_mod$variational(
  data = data_for_stan,
  seed = 123,
  iter = 40000
)

# Compute LOO
log_lik_vb <- sampled_vb$draws("log_lik")  
loo_vb <- loo(log_lik_vb)
print(loo_vb)

# Get the summary
vb_summary <- sampled_vb$summary()

# Convert the summary tibble to a data frame for easier indexing
vb_summary_df <- as.data.frame(vb_summary)

# Extract the 90% CIs for the mu_beta parameters (effects of the predictors)
ci_lower_mu_beta <- 
  vb_summary_df[grepl("mu_beta", vb_summary_df$variable), "q5"]
ci_upper_mu_beta <- 
  vb_summary_df[grepl("mu_beta", vb_summary_df$variable), "q95"]

# Create a data frame to neatly store these values
ci_df <- data.frame(
  Predictor = c("Outcome", "Stimulus", "RPE", "ZMoodPre", "ZControl", "Trial"),
  CI_Lower = ci_lower_mu_beta,
  CI_Upper = ci_upper_mu_beta
)

# Print the data frame
print(ci_df)



# Sample with MCMC
sampled <- sample_mod$sample(
  data = data_for_stan,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 300,
  iter_warmup = 100,
  refresh = 20
)

sampled$diagnostic_summary()
sampled$summary()

# Compute LOO
loo_m5 <- sampled$loo(cores = 4)
print(loo_m3)

summary_mod <- sampled$summary()
summary_mod$variable
summary(summary_mod$rhat)

# Get the summary
sampled_summary <- sampled$summary()

# Convert the summary tibble to a data frame for easier indexing
sampled_summary_df <- as.data.frame(sampled_summary)

# Extract the 95% CIs for the mu_beta parameters (effects of the predictors)
ci_lower_mu_beta <- sampled_summary_df[grepl("mu_beta", sampled_summary_df$variable), "q5"]
ci_upper_mu_beta <- sampled_summary_df[grepl("mu_beta", sampled_summary_df$variable), "q95"]

# Create a data frame to neatly store these values
ci_df <- data.frame(
  Predictor = c("Outcome", "Stimulus", "RPE", "ZMoodPre", "ZControl", "Trial"),
  CI_Lower = ci_lower_mu_beta,
  CI_Upper = ci_upper_mu_beta
)

# Print the data frame
print(ci_df)



# Plot -------------------------------------------------------------------------

# Extract the Samples
Y_pred_draws <- sampled_vb$draws(variables = "y_pred")
# Convert to an array
Y_pred_array <- as.array(Y_pred_draws)

# Compute the mode for each trial
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) mlv(as.vector(x)))
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) median(as.vector(x)))
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) mean(as.vector(x), trim=0.1))

Y_pred <- colMeans(Y_pred_array)

# Assume Y_observed is a vector of observed values
comparison <- data.frame(
  Y_observed = data_for_stan$happiness, 
  Y_pred = Y_pred
)

cor(
  comparison$Y_observed, comparison$Y_pred
)


all_data_df$happiness_hat <- Y_pred


# Calculate mean and standard error of the mean (SEM)
for_plot_df <- all_data_df |> 
  group_by(is_reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(happiness_hat, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(happiness_hat, na.rm = TRUE) / sqrt(n())
  ) |> 
  ungroup()
for_plot_df$is_reversal <- as.factor(for_plot_df$is_reversal)

# Create the plot
for_plot_df |> 
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = viridis(3)[1], linewidth = 2) +  # Line for h
  geom_ribbon(aes(ymin = h - h_sem, ymax = h + h_sem), alpha = 0.2) +  # Ribbon for h
  geom_line(aes(y=h_hat), color = viridis(3)[2], linewidth = 1.0) +  # Line for h_hat
  #geom_ribbon(aes(ymin = h_hat - h_hat_sem, ymax = h_hat + h_hat_sem), alpha = 0.2, fill = "red") +  # Ribbon for h_hat
  facet_wrap(~ is_reversal) +
  theme_default()  # Apply the bayesplot theme



## NOT BAD!!!! Tue Oct 31 21:34:29 2023


# Single session

prepare_stan_data_single_session <- function(
    df, target_user_id, target_ema_number) {
  
  # Filter data for a specific user_id and ema_number (session)
  filtered_df <- df %>%
    filter(user_id == target_user_id, ema_number == target_ema_number) %>%
    select(trial, happiness, RPE, outcome, stimulus, zmoodpre, zcontrol) %>%
    arrange(trial)
  
  # Number of observations in this session
  N <- nrow(filtered_df)
  
  # Maximum trial number in this session
  T <- max(filtered_df$trial)
  
  # Create the list to be passed to Stan
  stan_data <- list(
    N = N,
    T = T,
    outcome = as.vector(filtered_df$outcome),
    stimulus = as.vector(filtered_df$stimulus),
    RPE = as.vector(filtered_df$RPE),
    zmoodpre = as.vector(filtered_df$zmoodpre),
    zcontrol = as.vector(filtered_df$zcontrol),
    trial = as.vector(filtered_df$trial),
    happiness = as.vector(filtered_df$happiness)
  )
  
  return(stan_data)
}



# Identify unique user_ids
set.seed(1)
unique_user_ids <- unique(all_data_df$user_id)
# Randomly sample 2 user_ids
random_user_ids <- sample(unique_user_ids, 3)
random_user_ids
# Filter the original data frame to include only these 10 user_ids
filtered_df <- all_data_df %>% 
  filter(user_id %in% random_user_ids)

# data_for_stan <- prepare_stan_data(filtered_df)



stan_data_single_session <- prepare_stan_data_single_session(
  all_data_df, target_user_id = unique_user_ids[49], 
  target_ema_number = 3)


# Check with VI whether the model converges
sampled_vb <- sample_mod$variational(
  data = stan_data_single_session,
  seed = 42,
  iter = 40000
)


# Plot -------------------------------------------------------------------------

# Extract the Samples
Y_pred_draws <- sampled_vb$draws(variables = "y_pred")
# Convert to an array
Y_pred_array <- as.array(Y_pred_draws)

# Compute the mode for each trial
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) mlv(as.vector(x)))
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) median(as.vector(x)))
# Y_pred_mode <- apply(Y_pred_array, 3, function(x) mean(as.vector(x), trim=0.1))

Y_pred <- colMeans(Y_pred_array)

# Assume Y_observed is a vector of observed values
comparison <- data.frame(
  Y_observed = stan_data_single_session$happiness, 
  Y_pred = Y_pred
)

cor(comparison$Y_observed, comparison$Y_pred)
plot(1:30, comparison$Y_observed)
points(1:30, comparison$Y_pred, type="l")




# Calculate mean and standard error of the mean (SEM)
for_plot_df <- all_data_df |> 
  group_by(is_reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(happiness_hat, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(happiness_hat, na.rm = TRUE) / sqrt(n())
  ) |> 
  ungroup()
for_plot_df$is_reversal <- as.factor(for_plot_df$is_reversal)

# Create the plot
for_plot_df |> 
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = viridis(3)[1], linewidth = 2) +  # Line for h
  geom_ribbon(aes(ymin = h - h_sem, ymax = h + h_sem), alpha = 0.2) +  # Ribbon for h
  geom_line(aes(y=h_hat), color = viridis(3)[2], linewidth = 1.0) +  # Line for h_hat
  #geom_ribbon(aes(ymin = h_hat - h_hat_sem, ymax = h_hat + h_hat_sem), alpha = 0.2, fill = "red") +  # Ribbon for h_hat
  facet_wrap(~ is_reversal) +
  theme_default()  # Apply the bayesplot theme

