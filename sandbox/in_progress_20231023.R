
# Script for instant mood analysis
# Project: groundhog_day

# Get estimates of the momentary happiness model for each user_id and
# ema_number (nrow = 2131)
# The tar_load() command is necessary only for debugging.
tar_load(params_happiness_clean_df)
tar_load(prl_df)



tar_load(params_happiness_df)

# Calculate mean and standard error of the mean (SEM)
for_plot_df <- params_happiness_df |> 
  group_by(is_reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(predicted_happiness, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(predicted_happiness, na.rm = TRUE) / sqrt(n())
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



