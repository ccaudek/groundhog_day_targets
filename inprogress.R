#' Script for the momentary mood analysis
#' Project: groundhog day
#' 
#' To use this script, source _targets.R:
source("_targets.R")


# Momentary happiness ----------------------------------------------------------

# Plot momentary happiness as a function of trial, for reversal and no-reversal
# sessions. Adds the predictions of the momentary happiness model 

# snippets/plot_happiness_model_predicitons.R


# Piecewise Regression Analysis ------------------------------------------------

# snippets/piecewise_regression_analysis.R




stancode(fit3_rev)

#### Questo è per controllare che i valori predetti dal modello 2 siano 
#### coerenti con le medie riportate nella figura che motiva l'analisi statistica.

# Calculate the mean fitted value for each observation
average_fitted_values <- apply(foo, 2, mean)  # Apply the mean function across rows (MCMC samples)

# Add the average fitted values to the original data frame
bysubj_rev_df$average_fitted_happiness <- average_fitted_values

# Now, if you want to get the average predicted happiness for each level of 't', you can group by 't' and summarize.
average_fitted_per_t <- bysubj_rev_df %>%
  group_by(trial) %>%
  summarize(average_happiness = mean(average_fitted_happiness, na.rm = TRUE))

average_fitted_per_t |> 
  ggplot(aes(x = trial, y = average_happiness)) +
  geom_line() +
  ggtitle("Average Predicted Happiness as a Function of t")

# ------------------------------------------------------------------------------

# The difference mood_post - mood_pre is accounted for by the regression 
# toward the mean.

temp <- params_happiness_clean_df
temp$mdif <- temp$mood_post - temp$mood_pre

bysubj_mood_df <- temp |> 
  group_by(user_id) |> 
  summarize(
    mood_dif = mean(mdif),
    m_pre = mean(mood_pre)
  ) |> 
  ungroup()

bysubj_mood_df$md <- as.vector(scale(bysubj_mood_df$mood_dif))
bysubj_mood_df$mp <- as.vector(scale(bysubj_mood_df$m_pre))

  
fm <- lm(
  md ~ mp,
  data = bysubj_mood_df
)
summary(fm)

mod1 <- brm(
  md ~ se * mp + (se * mp | user_id),
  data = bysubj_sess_mood_df,
  backend = "cmdstanr",
  cores = 8, 
  chains = 4,
  iter = 5000,
  prior = c(
    set_prior("normal(0, 1)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 1)", class = "sd"), # Weakly informative prior for random effects sd
    set_prior("lkj(2)", class = "cor") # Weakly informative prior for random effects correlations
  ),
  threads = threading(2),
)

summary(mod1)

conditional_effects(mod1, "mp:se")
performance::r2_bayes(mod1)


# ------------------------------------------------------------------------------

# This works, but requires more iterations.
fit_rev <- brm(
  happiness ~ epoch + t + (epoch + t | user_id / ema_number),
  data = rev_df,
  family = asym_laplace(),
  backend = "cmdstanr",
  prior = c(
    set_prior("normal(0, 0.5)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 1)", class = "sd", group = "user_id"), # Weakly informative prior for random effects sd
    set_prior("lkj(2)", class = "cor", group = "user_id") # Weakly informative prior for random effects correlations
  ),
  cores = 4,
  chains = 4,
  threads = threading(4)
)
pp_check(fit_rev)

summary(fit_rev)
conditional_effects(fit_rev, "t")
bayes_R2(fit_rev)

performance::r2_bayes(fit_rev)

loo_fit_rev <- loo(fit_rev)








mod <- cmdstan_model("R/stan_code/M11.stan", compile = FALSE)
mod$format(canonicalize = TRUE)

# Assuming each subject has a unique user_id and each session is uniquely identified by a combination of user_id and ema_number
# Create a new column to uniquely identify sessions
rev_df <- rev_df %>%
  mutate(session_id = interaction(user_id, ema_number, drop = TRUE))

# Create a vector to indicate which subject each observation belongs to
subject <- as.integer(factor(rev_df$user_id))

# Create a vector to indicate which session each observation belongs to
session <- as.integer(factor(rev_df$ema_number))

# Get the total number of observations, subjects, and trials per session
N <- nrow(rev_df)
S <- length(unique(rev_df$user_id))
T <- 30  # since you mentioned each session has 30 trials

# Get the number of sessions per subject
n_sessions <- table(rev_df$user_id)

# Prepare the data list for Stan
stan_data <- list(
  N = N,
  S = S,
  T = T,
  n_sessions = n_sessions,
  x = rev_df$trial,
  y = rev_df$happiness,  # replace with the name of your response variable
  subject = subject,
  session = session
)

# Now you can use stan_data as the data input for your Stan model


# Compile M9_piecewise.stan
sample_mod <- cmdstan_model(
  "R/stan_code/M11.stan", 
  stanc_options = list("O1"),
  force_recompile = TRUE
)

# Use Variational Inference
sampled_vb <- sample_mod$variational(
  data = stan_data,
  seed = 123,
  iter = 40000
)
sampled_vb$summary()

# Sample with MCMC
sampled <- sample_mod$sample(
  data = stan_data,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 700,
  iter_warmup = 200,
  refresh = 20
)
sampled$summary()

beta1 <- sampled$summary("beta1")
beta2 <- sampled$summary("beta2")
mu_beta1 <- sampled$summary("mu_beta1")
mu_beta2 <- sampled$summary("mu_beta2")

draws <- sampled$draws(format = "df")

df <- data.frame(
  beta = c(draws$mu_beta1, draws$mu_beta2),
  segment = factor(rep(c("Segment 1", "Segment 2"), each = length(draws$mu_beta1)))
)

ggplot(df, aes(x = segment, y = beta, fill = segment)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(
    title = "Group Effects for mu_beta1 and mu_beta2",
    y = "Parameter Estimate",
    x = "Segment"
  )









for_plot_df <- bysubj_rev_df |> 
  group_by(trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1)
  )

plot(for_plot_df$trial, for_plot_df$h, type = 'l')



fit_norev <- brm(
  happiness ~ 1 + trial + (1 + trial | user_id / ema_number),
  data = norev_df,
  family = gaussian(),
  backend = "cmdstanr",
  chains = 2,
  iter = 500,
  prior = c(
    set_prior("normal(0, 5)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 2)", class = "sd") # Weakly informative prior for random effects
  ),
  # decomp = "QR",
  # normalize = FALSE,
  threads = threading(4)
)
pp_check(fit_norev) #+ xlim(-3, 3)

summary(fit_norev)
marginal_effects(fit_norev, "trial")

performance::r2_bayes(fit_norev)





bysubj_rev_df <- rev_df |> 
  group_by(user_id, trial) |> 
  summarize(
    happiness = mean(happiness, trim = 0.1),
  ) |> 
  ungroup()

for_plot_df <- bysubj_rev_df |> 
  group_by(trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1)
  )

plot(for_plot_df$trial, for_plot_df$h, type = 'l')

bysubj_rev_df$epoch <- ifelse(
  bysubj_rev_df$trial < 16, 0, 1
) |> 
  as.factor()

fit_rev <- brm(
  happiness ~ 1 + epoch * trial + (1 + epoch * trial | user_id),
  data = bysubj_rev_df,
  family = student(),
  backend = "cmdstanr",
  chains = 2,
  iter = 1000,
  prior = c(
    set_prior("normal(0, 5)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 5)", class = "sd") # Weakly informative prior for random effects
  ),
  threads = threading(4)
)
pp_check(fit_rev) + xlim(-3, 3)

summary(fit_rev)
marginal_effects(fit_rev, "trial")
marginal_effects(fit_rev, "epoch")
marginal_effects(fit_rev, "trial:epoch")


r2_bayes(fit_rev)


# Piecewise regression with Stan


# Center the trial variable around the change point
bysubj_rev_df$centered_trial <- bysubj_rev_df$trial - 16

# Assuming 'centered_trial' is centered around 0, we create the new variables 
# for the slopes
bysubj_rev_df$trial_pre = 
  ifelse(bysubj_rev_df$centered_trial < 0, bysubj_rev_df$centered_trial, 0)
bysubj_rev_df$trial_post = 
  ifelse(bysubj_rev_df$centered_trial >= 0, bysubj_rev_df$centered_trial, 0)















# Create an indicator variable for trials before and after the change point
bysubj_rev_df$pre_change <- ifelse(bysubj_rev_df$trial < 16, 1, 0)
bysubj_rev_df$post_change <- ifelse(bysubj_rev_df$trial >= 16, 1, 0)

# Create the design matrix with centered trial
# Since we're dealing with a piecewise regression, we don't need interaction with 'epoch'
X_pre <- model.matrix(~ centered_trial * pre_change, bysubj_rev_df)
X_post <- model.matrix(~ centered_trial * post_change, bysubj_rev_df)

# Prepare the data for Stan
stan_data <- list(
  N = nrow(bysubj_rev_df),
  Y = bysubj_rev_df$happiness,
  trial = bysubj_rev_df$centered_trial,
  N_subjects = length(unique(bysubj_rev_df$user_id)),
  subject = as.integer(factor(bysubj_rev_df$user_id))  # Convert user_id to consecutive integers
)


# Compile M9_piecewise.stan
sample_mod <- cmdstan_model(
  "R/stan_code/M10.stan", 
  stanc_options = list("O1"),
  force_recompile = TRUE
)

# Use Variational Inference
sampled_vb <- sample_mod$variational(
  data = stan_data,
  seed = 123,
  iter = 40000
)

# Sample with MCMC
sampled <- sample_mod$sample(
  data = stan_data,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 500,
  iter_warmup = 100,
  refresh = 20
)

vb_summary <- sampled_vb$summary()
vb_summary$variable

draws <- sampled_vb$draws(format = "df")


# Identify all y_pred columns
y_pred_columns <- grep("y_pred", names(draws), value = TRUE)

# Pivot the data from wide to long format
draws_long <- draws %>%
  pivot_longer(cols = y_pred_columns, names_to = "trial", values_to = "y_pred") %>%
  mutate(trial = readr::parse_number(trial)) # Extract the numeric part of the trial identifier



predicted_Y <- draws$

# Determine the number of draws per observation
num_draws_per_obs <- nrow(draws) / nrow(stan_data)

# Repeat the stan_data for each draw
stan_data_long <- stan_data[rep(1:nrow(stan_data), each = num_draws_per_obs),]

# Now combine the repeated stan_data with the draws
combined_data <- cbind(stan_data_long, draws)




# Create a new column for the combined effect
draws$combined_effect <- draws$`b[2]` + draws$`b[3]`
hist(draws$combined_effect, breaks=50, main="Posterior of Combined Effect of Trial after Change Point")

hist(draws$intercept)


draws <- sampled$draws(format = "df")

# Create a new column for the combined effect
draws$combined_effect <- draws$`b[2]` + draws$`b[3]`
hist(draws$combined_effect, breaks=50, main="Posterior of Combined Effect of Trial after Change Point")








# Define the known change point
known_omega <- 16

bform <- bf(
  y ~ b0 + b1 * (trial - known_omega) * step(known_omega - trial) + 
    b2 * (trial - known_omega) * step(trial - known_omega),
  b0 + b1 + b2 ~ 1 + (1 | user_id)
)

df <- data.frame(
  y = rnorm(330),
  age = rep(0:10, 30),
  person = rep(1:30, each = 11)
)

bprior <- prior(normal(0, 3), nlpar = "b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2")

make_stancode(bform, data = df, prior = bprior)







# Relative importance of the predictors of momentary mood

mydat <- params_happiness_clean_df 
mydat$happiness <- ifelse(
  mydat$happiness < -5 | mydat$happiness > 5, NA, mydat$happiness)


foo <- mydat |> 
  dplyr::select(c(w_outcome, w_stimulus, w_rpe, w_moodpre, w_control, w_trial))

cor(foo)


norev_df <- params_happiness_clean_df |> 
  dplyr::filter(is_reversal == 0) |> 
  group_by(user_id, trial) |> 
  summarize(
    happiness = mean(happiness), 
    w_outcome = mean(w_outcome),
    w_rpe = mean(w_rpe),
    w_moodpre = mean(w_moodpre),
    w_control = mean(w_control),
    w_trial = mean(w_trial)
  ) |> 
  ungroup()

for_plot_df <- norev_df |> 
  group_by(trial) |> 
  summarize(
    h = mean(happiness),
    w_trial = mean(w_trial)
  )

plot(for_plot_df$w_trial, for_plot_df$h, type = 'l')

fit_norev <- brm(
  h ~ trial + 
    w_outcome + w_stimulus + w_rpe + w_moodpre + w_control +
    (is_reversal * w_trial + w_outcome + w_stimulus + w_rpe + w_moodpre + 
       w_control | user_id / ema_number),
  data = mydat,
  family = asym_laplace(),
  prior = c(
    set_prior("normal(0, 5)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 2)", class = "sd") # Weakly informative prior for random effects
  ),
  algorithm = "meanfield"
)
pp_check(fit2)





fit2 <- brm(
  happiness ~ is_reversal * w_trial + 
    w_outcome + w_stimulus + w_rpe + w_moodpre + w_control +
    (is_reversal * w_trial + w_outcome + w_stimulus + w_rpe + w_moodpre + 
       w_control | user_id / ema_number),
  data = mydat,
  family = asym_laplace(),
  prior = c(
    set_prior("normal(0, 5)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 2)", class = "sd") # Weakly informative prior for random effects
  ),
  algorithm = "meanfield"
)
pp_check(fit2)

summary(fit2)

bayes_R2(fit2)
#     Estimate   Est.Error      Q2.5     Q97.5
# R2 0.4555724 0.009869782 0.4365242 0.4744727

conditional_effects(fit2, "w_trial")
loo_fit2 <- loo(fit2)


fit3 <- brm(
  happiness ~ is_reversal + 
    w_outcome + w_stimulus + w_rpe + w_moodpre + w_control +
    (is_reversal + w_outcome + w_stimulus + w_rpe + w_moodpre + 
       w_control | user_id / ema_number),
  data = mydat,
  family = asym_laplace(),
  prior = c(
    set_prior("normal(0, 5)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 2)", class = "sd") # Weakly informative prior for random effects
  ),
  algorithm = "meanfield"
)
bayes_R2(fit3)

loo_fit3 <- loo(fit3)

loo_compare(loo_fit2, loo_fit3)


# Using Stan

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
library(rstan)
library(loo)
color_scheme_set("brightblue")

cmdstan_version()


# Compile M8_picewise.stan
sample_mod <- cmdstan_model(
  "R/stan_code/M8_picewise.stan", 
  stanc_options = list("O1"),
  force_recompile = TRUE
)

# Create list for Stan data
data_for_stan <- list(
  N = nrow(params_happiness_clean_df),
  happiness = params_happiness_clean_df$happiness,
  X = as.matrix(params_happiness_clean_df[, c('w_outcome', 'w_stimulus', 'w_rpe', 'w_moodpre', 'w_control', 'w_trial')]),
  J = length(unique(params_happiness_clean_df$user_id)),
  subject = as.integer(factor(params_happiness_clean_df$user_id)),
  is_reversal = as.integer(params_happiness_clean_df$is_reversal),
  trial = params_happiness_clean_df$trial
)

str(stan_data)


# Function to generate initial values
init_function <- function() {
  return(list(
    alpha = runif(1, 0.1, 0.9),  # More constrained
    beta = runif(6, 1.1, 1.9),  # More constrained and assuming beta is a vector of length 6
    beta_reversal = runif(1, 0.1, 0.9),
    beta_piecewise1 = runif(1, 0.1, 0.9),
    beta_piecewise2 = runif(1, 0.1, 0.9),
    sigma = runif(1, 0.1, 0.9),
    sigma_subject = runif(1, 0.1, 0.9),
    z_subject_raw = rnorm(215, 0, 1),  # Centered around 0 with std dev 1
    kappa = runif(1, 0.1, 0.9)
  ))
}


# Use Variational Inference
sampled_vb <- sample_mod$variational(
  data = data_for_stan,
  # init = init_function,
  seed = 123,
  iter = 40000
)

draws_df <- as.data.frame(sampled_vb$draws())


# Extract alpha and beta samples
alpha_samples <- draws_df$alpha
beta_samples <- draws_df[, grep("beta", names(draws_df))]

# Convert beta_samples to a matrix
beta_samples_matrix <- as.matrix(beta_samples)

# Compute the posterior predictive mean for each observation
posterior_predictive_means <- matrix(0, nrow = nrow(X_matrix), ncol = length(alpha_samples))

for (s in 1:length(alpha_samples)) {
  posterior_predictive_means[, s] <- alpha_samples[s] + X_matrix %*% beta_samples_matrix[s,]
}

# Compute the mean and 95% CI for the posterior predictive distribution
post_pred_mean <- rowMeans(posterior_predictive_means)
post_pred_lower <- apply(posterior_predictive_means, 1, function(x) quantile(x, 0.025))
post_pred_upper <- apply(posterior_predictive_means, 1, function(x) quantile(x, 0.975))

# Plot the observed data and the posterior predictive checks
plot(data_for_stan$happiness, type = 'l', col = 'black', lty = 1, ylim = c(min(post_pred_lower), max(post_pred_upper)), xlab = "Index", ylab = "Happiness")
lines(post_pred_mean, col = 'blue', lty = 2)
lines(post_pred_lower, col = 'red', lty = 2)
lines(post_pred_upper, col = 'red', lty = 2)
legend("topright", legend = c("Observed", "Predicted", "95% CI"), col = c("black", "blue", "red"), lty = c(1, 2, 2))


data_for_stan$predicted_happiness <- post_pred_mean


# Aggregate observed data
agg_data_observed <- params_happiness_clean_df %>%
  group_by(is_reversal, trial) %>%
  summarise(mean_happiness = mean(happiness)) |> 
  ungroup()

params_happiness_clean_df$predicted_happiness <- post_pred_mean

agg_data_predicted <- params_happiness_clean_df %>% 
  group_by(is_reversal, trial) %>%
  summarise(mean_predicted_happiness = mean(predicted_happiness)) |> 
  ungroup()

ggplot() +
  #geom_line(data = agg_data_observed, aes(x = trial, y = mean_happiness), color = "black", alpha = 0.5) +
  geom_line(data = agg_data_predicted, aes(x = trial, y = mean_predicted_happiness), color = "red", alpha = 0.5) +
  facet_wrap(~ is_reversal, ncol = 2) +
  labs(title = "Model Fit by Session Type",
       x = "Trial",
       y = "Mean Happiness") +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_minimal()

#####

# Step 1: Create a grid for 'trial'
trial_grid = seq(min(params_happiness_clean_df$trial), max(params_happiness_clean_df$trial), length.out = 100)

# Step 2 & 3: Compute Predicted Happiness and Average Them
pdp_values = numeric(length(trial_grid))

for (i in seq_along(trial_grid)) {
  trial_value = trial_grid[i]
  
  # Create a modified data frame where all 'trial' values are set to the current grid value
  modified_data = params_happiness_clean_df
  modified_data$trial = trial_value
  
  # Extract the weights for this grid value (assuming weights are constant across trials)
  params = c(
    modified_data$w0[1], modified_data$w_outcome[1], modified_data$w_stimulus[1], 
    modified_data$w_rpe[1], modified_data$w_moodpre[1], modified_data$w_control[1], 
    modified_data$w_trial[1], modified_data$w_gamma[1]
  )
  
  # Compute the predicted happiness using your original function
  predicted_happiness = compute_predicted_happiness(params, modified_data)
  
  # Average these predictions to get the partial dependence for this 'trial' value
  pdp_values[i] = mean(predicted_happiness)
}

# Step 4: Plot the Partial Dependence Plot
plot(trial_grid, pdp_values, type = 'l', xlab = 'Trial', ylab = 'Average Predicted Happiness')

####

library(relaimpo)

rw <- calc.relimp(
  happiness ~ w_outcome + w_stimulus + w_rpe + w_moodpre + w_control + 
    w_trial, 
  data = params_happiness_clean_df, 
  type = c("lmg"),
  rela=TRUE
  )
print(rw)

library(randomForest)

# Selecting predictors and the dependent variable
predictors <- params_happiness_clean_df[, c("w_outcome", "w_stimulus", "w_rpe", "w_moodpre", "w_control", "w_trial")]
response <- params_happiness_clean_df$happiness

# Train Random Forest model
rf_model <- randomForest(x=predictors, y=response, ntree=500)

# Show model summary
print(rf_model)

# Extract variable importance
importance(rf_model)





# eof ----














# Ungroup (otherwise there will be an error!) and change name 
# to the data frame. 
params <- params_happiness_clean_df |> 
  ungroup()

# Ci sono 22 soggetti con alpha = 0. Andrebbero rimossi.

# Initialize an empty data frame to store the results
final_df <- data.frame()

# prl_df is made available in the workflow. 
# Compute `zim`: the perceived mood in each trial of the APRL task, which is
# standardized for each participant.
dz_clean <- prl_df |> 
  group_by(user_id) |> 
  mutate(
    zim = as.vector(scale(instant_mood))
  )

# Get all unique user_ids.
unique_users <- unique(dz_clean$user_id)

# Loop through each unique user_id.
for (id in unique_users) {
  # Select one subject 
  onesubj_data <- dz_clean %>% dplyr::filter(user_id == id)
  
  # Get alpha for all ema_number sessions and the current user_id
  best_alpha = get_alpha_all_sessions(onesubj_data)
  
  n_ema_episodes <- length(unique(onesubj_data$ema_number))
  
  # Loop through each ema_number for the current user_id.
  for (i in seq_len(n_ema_episodes)) {
    ema_session <- onesubj_data %>% dplyr::filter(ema_number == i)
    
    # Required information for a single session of a subject.
    df <- data.frame(
      trial = ema_session$trial,
      stimulus = ifelse(
        ema_session$is_target_chosen == 0, -1, ema_session$is_target_chosen
      ), # 1 for stimulus A, -1 for stimulus B
      reversal = ifelse(
        ema_session$is_reversal == "yes",
        c(rep(0, 15), 1, rep(0, 14)),
        rep(0, 30)
      ),
      outcome = ifelse(ema_session$feedback == 0, -1, ema_session$feedback),
      delta_p = ifelse(
        ema_session$is_reversal == "yes",
        c(rep(0, 15), 0.6, rep(0, 14)),
        rep(0, 30)
      ),
      happiness = ema_session$zim, # standardized by user_id
      zmoodpre = ema_session$zmoodpre,
      zcontrol = ema_session$zcontrol
    )
    
    # Add the RPE column to the df DataFrame
    df = add_rpe(df, best_alpha)
    
    subj_code <- unique(ema_session$user_id)
    
    subj_session_params <- params |> 
      dplyr::filter(user_id == subj_code & ema_number == i) |> 
      dplyr::select("w0", "w1", "w2", "w3", "w4", "w5", "w6", "gamma")
    
    happiness_hat <- compute_predicted_happiness(
      subj_session_params, df
    ) |> unlist()
    
    # After calculating happiness_hat, create a temporary data frame
    temp_df <- data.frame(
      user_id = id,
      ema_session = i,
      trial = df$trial,
      reversal = ema_session$is_reversal,
      stimulus = df$stimulus,
      outcome = df$outcome,
      rpe = df$RPE,
      best_alpha = best_alpha,
      happiness = df$happiness,
      happiness_hat = happiness_hat
    )
    
    # Append this data to the final data frame
    final_df <- rbind(final_df, temp_df)
  }
}


glimpse(final_df)

hist(final_df$happiness)

hist(final_df$happiness_hat)

# Whether alpha is below 0.05 or above does not change the fit of
# happiness_hat to happiness.
good_alpha_df <- final_df[final_df$best_alpha > 0.05, ]
bad_alpha_df <- final_df[final_df$best_alpha <= 0.05, ]


# Calculate mean and standard error of the mean (SEM)
for_plot_df <- final_df |> 
  # dplyr::filter(ema_session == 12) |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(happiness_hat, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(happiness_hat, na.rm = TRUE) / sqrt(n())
  ) |> 
  ungroup()
for_plot_df$reversal <- as.factor(for_plot_df$reversal)

# Create the plot
for_plot_df |> 
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = viridis(3)[1]) +  # Line for h
  geom_ribbon(aes(ymin = h - h_sem, ymax = h + h_sem), alpha = 0.2) +  # Ribbon for h
  geom_line(aes(y=h_hat), color = viridis(3)[2]) +  # Line for h_hat
  geom_ribbon(aes(ymin = h_hat - h_hat_sem, ymax = h_hat + h_hat_sem), alpha = 0.2, fill = "red") +  # Ribbon for h_hat
  facet_wrap(~ reversal) +
  theme_default()  # Apply the bayesplot theme


# Plot for a single user_id
for_plot_df <- final_df |> 
  dplyr::filter(user_id == 3881466599) |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(happiness_hat, trim = 0.1, na.rm = TRUE)
  ) |> 
  ungroup()
for_plot_df$reversal <- as.factor(for_plot_df$reversal)

# Create the plot
for_plot_df %>%
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = viridis(3)[1]) +  # Line for h
  geom_line(aes(y=h_hat), color = viridis(3)[2]) +  # Line for h_hat
  facet_wrap(~ reversal) +
  theme_default()  # Apply the bayesplot theme


# Correlation empirical vs predicted happiness
cor.test(
  final_df$happiness, final_df$happiness_hat, 
  na.action=na.omit, method = "pearson"
)










fm <- brm(
  mood_dif ~ mood_pre + environment + (w0 + w1 + w2 + w3 + w4 + w5 + w6) +
    (mood_pre + environment | user_id / ema_number),
  data = params_happiness_clean_df,
  family = student(),
  algorithm = "meanfield"
  #backend = "cmdstanr"
)
pp_check(fm) + xlim(-80, 80)
bayes_R2(fm)
performance::r2_bayes(fm)
summary(fm)
conditional_effects(fm, "w0")
conditional_effects(fm, "environment")
conditional_effects(fm, "mood_pre")


performance::r2_bayes(fm)

temp <- prl_df

temp1 <- temp[temp$is_reversal == "no" & temp$ema_number < 13, ]
temp1$ses <- as.vector(scale(temp1$ema_number))
temp1$tr <- as.vector(scale(temp1$trial))
temp1$im <- as.vector(scale(temp1$instant_mood))


bysubj_params_df <- params_happiness_clean_df |>
  dplyr::select(
    "w0", "w1", "w2", "w3", "w4",
    "w5", "w6", "gamma", "is_reversal", "ema_number",
    "user_id", "alpha", "environment"
  ) 
bysubj_params_df$is_reversal <- ifelse(
  bysubj_params_df$is_reversal == 1, "yes", "no"
)

mydat <- left_join(
  prl_df, bysubj_params_df, by = c("user_id", "ema_number", "is_reversal")
)


stable_df <- mydat |> 
  dplyr::filter(is_reversal == "no")

volatile_df <- mydat |> 
  dplyr::filter(is_reversal == "yes")

volatile_df$change_point <- ifelse(volatile_df$trial > 15, 1, 0)

fit0 <- brm(
  instant_mood ~ 1 * change_point + trial * change_point +
    (1 + trial | user_id / ema_number), 
  data = volatile_df,
  family = student(), #cumulative("logit"),
  #algorithm = "meanfield",
  backend = "cmdstanr",
  cores = 4
  #iter = 4000,
)

pp_check(fit0) # + xlim(-4, 4)
bayes_R2(fit0)
summary(fit0)
conditional_effects(fit0, "trial")


fit1 <- brm(
  instant_mood ~ 1 + trial + 
    (1 + trial | user_id / ema_number), 
  data = stable_df,
  family = cumulative("logit"),
  algorithm = "meanfield",
  # backend = "cmdstanr",
  cores = 4
)

pp_check(fit1) # + xlim(-4, 4)
bayes_R2(fit1)
summary(fit1)
conditional_effects(fit1, "trial")
conditional_effects(fit1, "tr")


fit2 <- brm(
  instant_mood ~ 1 + 
    w0 + w1 + w2 + w3 + w4 + w5 + w6 + 
    (1 + trial | user_id / ema_number / trial), 
  data = stable_df,
  family = cumulative("logit"),
  algorithm = "meanfield",
  # backend = "cmdstanr",
  cores = 4
)

pp_check(fit2) # + xlim(-4, 4)
bayes_R2(fit2)
summary(fit2)
conditional_effects(fit2, "w0")
conditional_effects(fit2, "w1")
conditional_effects(fit2, "w3")
conditional_effects(fit2, "w4")
conditional_effects(fit2, "w5")
conditional_effects(fit2, "w6")


stable_emanum_df <- stable_df |> 
  group_by(user_id, ema_number) |> 
  summarize(
    instant_mood = mean(instant_mood),
    w0 = mean(w0), 
    w1 = mean(w1), 
    w2 = mean(w2), 
    w3 = mean(w3), 
    w4 = mean(w4), 
    w5 = mean(w5), 
    w6 = mean(w6),
    gamma = mean(gamma)
  ) |> 
  ungroup()

# Only result: the weight given to the outcomes increaes with sessions.

fit4 <- brm(
  w0 ~ 1 + ema_number + 
    (1 + ema_number | user_id), 
  data = mydat[mydat$is_reversal == "no", ],
  family = student(),
  algorithm = "meanfield",
  # backend = "cmdstanr",
  cores = 4
)
pp_check(fit4) + xlim(-2, 2)
bayes_R2(fit4)
summary(fit4)
conditional_effects(fit4, "ema_number")


delta_t <-
  # extracting posterior samples from bmod5
  posterior_samples(fit1, pars = c("^b_", "sd_", "sigma") ) %>% # taking the square of each variance component
  mutate_at(.vars = 6:9, .funs = funs(.^2) ) %>%
  # dividing the slope estimate by the square root of the sum of # all variance components
  mutate(delta = b_tr / sqrt(rowSums(.[6:9]) ) )

mean(delta_t$delta)








temp1 <- prl_df

temp1$ses <- as.vector(scale(temp1$ema_number))
temp1$tr <- as.vector(scale(temp1$trial))
temp1$im <- as.vector(scale(temp1$instant_mood))

temp1$tr_segment1 <- ifelse(temp1$is_reversal == 1 & temp1$tr <= 15, temp1$tr, 0)
temp1$tr_segment2 <- ifelse(temp1$is_reversal == 1 & temp1$tr > 15, temp1$tr - 15, 0)



fit1 <- brm(
  instant_mood ~ 
    1 +
    tr * (1 - is_reversal) +  # Effect of tr when is_reversal is 0
    tr_segment1 * is_reversal +  # Effect of tr_segment1 when is_reversal is 1
    tr_segment2 * is_reversal +  # Effect of tr_segment2 when is_reversal is 1
    (1 + tr | user_id / ses),
  data = temp1,
  family = cumulative("logit"),
  algorithm = "meanfield",
  cores = 4,
  prior = c(
    prior(normal(0, 4), class = "Intercept"),
    prior(normal(0, 4), class = "b")
  )
)






# Extract posterior samples
post_samples <- posterior_samples(fit1, pars = "b_tr")


# Calculate the mean of the posterior samples for 'tr'
mean_tr <- mean(post_samples$b_tr)

# Calculate the pooled standard deviation of the outcome variable
# You should calculate this on your original data (temp1 in your case)
pooled_sd <- sd(temp1$instant_mood)

# Calculate Cohen's d
cohen_d <- mean_tr / pooled_sd
cohen_d

#### QUEST

scs_scores <- rio::import(
  here::here(
    "data", "prep", "quest_scales", "scs_scores.csv"
  )
)
scs_scores$user_id <- as.numeric(scs_scores$user_id)

mydat <- left_join(params_happiness_clean_df, scs_scores, by = "user_id")

mydat <- mydat[!is.na(mydat$scs_total_score), ]
length(unique(mydat$user_id))

mydat$scs <- as.vector(scale(mydat$scs_total_score))
mydat$mpre <- as.vector(scale(mydat$mood_pre))
mydat$mpost <- as.vector(scale(mydat$mood_post))




######################################
# Outliers detection

params_values <- params_happiness_df |> 
  dplyr::select(w0, w1, w2, w3, w4, w5, gamma)

# Finding the center point 
params_center <- colMeans(params_values)

# Finding the covariance matrix
params_cov <- cov(params_values)

# Calculate Mahalanobis distance and add to data
params_values$mdist <- mahalanobis(
  x = params_values,
  center = params_center,
  cov = params_cov
)

# Cutoff values from chi-square distribution to identify outliers
cutoff <- qchisq(p = 0.99, df = ncol(params_values[, 1:7]))

# Add outliers based on Mahalanobis distance to data
params_values <- params_values %>%
  mutate(is_outlier = ifelse(mdist > cutoff, 1, 0))



table(params_values$is_outlier)

foo <- params_values[params_values$is_outlier == 0, ]





fit <- brm(
  happiness ~ happiness_hat + (1 + happiness_hat | user_id) + 
    (1 + happiness_hat | user_id:ema_session),
  data = final2_df,
  family = student(),  # or other appropriate family
  algorithm = "meanfield"
)

bayes_R2(fit)


fit <- brm(
  zim ~ zh  + 
    (zh | user_id / ema_session/ trial),
  data = temp,
  family = asym_laplace(),  
  algorithm = "meanfield",
  iter = 40000
)



# Correlation between momentary earning and momentary happiness
temp <- final2_df
temp <- temp |> 
  group_by(user_id, ema_session) |> 
  mutate(
    earning = cumsum(outcome)
  )
cor(temp$happiness, temp$earning, use="pairwise.complete.obs")


temp$zh <- scale(temp$h_hat)
temp$ze <- scale(temp$earning)

fit <- brm(
  happiness ~ zh + ze + 
    (zh + ze | user_id / ema_session),
  data = temp,
  family = asym_laplace(),  
  algorithm = "meanfield",
  iter = 40000
)

summary(fit)

posterior_samples_fit <- posterior_samples(fit)

h_hat_samples <- posterior_samples_fit$b_h_hat
earning_samples <- posterior_samples_fit$b_earning

diffs <- h_hat_samples - earning_samples

proportion_greater <- mean(diffs > 0)
proportion_smaller <- mean(diffs < 0)

# If proportion_greater is close to 1, then you have strong evidence that 
# h_hat is more strongly associated with happiness than earning is.

n_bootstraps <- 5000
bootstrap_results <- numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  bootstrap_sample <- sample(diffs, size = length(diffs), replace = TRUE)
  bootstrap_results[i] <- mean(bootstrap_sample > 0)
}

lower_bound <- quantile(bootstrap_results, 0.025)
upper_bound <- quantile(bootstrap_results, 0.975)

c(lower_bound, upper_bound)


mod_happiness <- brm(
  h ~ h_hat + reversal * (trial + ema_session) +
    (happiness_hat + reversal * (trial + ema_session) | user_id),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = final_df
)

pp_check(mod_happiness) + xlim(-10, 10)
bayes_R2(mod_happiness)
summary(mod_happiness)
marginal_effects(mod_happiness, "h_hat")

predicted_values <- posterior_predict(mod_happiness, newdata = final_df)

# The output will be a matrix where each row corresponds to an observation in final_df
# and each column corresponds to a posterior draw. You may want to summarize this 
# into a single predicted value per observation, e.g., by taking the mean across columns.

# Take the mean across the columns to get a single predicted value per observation
mean_predicted_values <- colMeans(predicted_values)


final_df |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness)
  )

foo <- final_df |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1),
    h_hat = mean(happiness_hat, trim = 0.1)
  )
foo$reversal <- as.factor(foo$reversal)

foo |> 
  ggplot(aes(x=trial, y=h, color = reversal)) +
  geom_line()


foo <- final_df |> 
  group_by(reversal, ema_session) |> 
  summarize(
    h = mean(happiness, trim = 0.1),
    h_hat = mean(happiness_hat, trim = 0.1)
  )
foo$reversal <- as.factor(foo$reversal)

foo |> 
  ggplot(aes(x=ema_session, y=h, color = reversal)) +
  geom_line()



