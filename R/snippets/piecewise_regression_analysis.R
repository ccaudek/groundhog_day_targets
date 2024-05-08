# Piecewise Regression Analysis ------------------------------------------------
#
# Date: Thu Nov 9 2023
# This script examines the increase of momentary happiness within a session of 
# 30 trials
# using a piecewise regression model with separate intercepts and slopes for 
# two epochs.
# The following implementation of the regression model with separate intercepts
# and slopes for the two epochs is working. Probably requires 10000 iterations.
# However, model comparison does not work.
# ELPD_DIFF does not show an advantage for the model with 'trial' as predictor
# of momentary happiness. 


# Load cleaned parameters for happiness
tar_load(params_happiness_clean_df)

# Split data by reversal status
# No reversal sessions
norev_df <- params_happiness_clean_df |> 
  dplyr::filter(is_reversal == 0) 

# Sessions with reversal
rev_df <- params_happiness_clean_df |> 
  dplyr::filter(is_reversal == 1) 

# Generate a data frame by collapsing over sessions, considering only reversal 
# sessions and summarizing happiness scores by user and trial, trimming 10% of 
# data
bysubj_rev_df <- rev_df |> 
  group_by(user_id, trial) |> 
  summarize(
    happiness = mean(happiness, trim = 0.1),
  ) |> 
  ungroup()

# Normalize trial numbers and create epoch variable
bysubj_rev_df$t <- as.vector(scale(bysubj_rev_df$trial))
bysubj_rev_df$epoch <- ifelse(bysubj_rev_df$trial < 16, 0, 1)
bysubj_rev_df$epoch <- factor(bysubj_rev_df$epoch)

# Plot mean happiness by trial to check for correctness
for_plot_df <- bysubj_rev_df |> 
  group_by(trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1)
  )
plot(for_plot_df$trial, for_plot_df$h, type = 'l')

# Frequentist check for the effect of time (t) using linear mixed-effects model
fit_lmer <- lmer(
  happiness ~ epoch * t + (epoch * t | user_id),
  data = bysubj_rev_df
)
summary(fit_lmer)

# To extract the random effects:
re <- ranef(fit_lmer)$user_id
# If you want to combine them with the fixed effects to get the total effect for each user:
fe <- fixef(fit_lmer)  # Extract the fixed effects
# Total slope for t for each user_id will be the fixed effect for t plus the random effect for t:
total_slopes_t <- fe["t"] + re[, "t"]
hist(total_slopes_t)
mean(total_slopes_t > 0)
# [1] 0.7627907

# Considering potential outliers, Student's t-distribution is used for the 
# response
plot(density(bysubj_rev_df$happiness))

# Fit the first model with the 'brm' function from the 'brms' package
# Includes both fixed effects and varying intercepts and slopes by user_id
# using a Student's t-distribution for the response
fit_rev <- brm(
  bf(happiness ~ 1 + epoch * t + 
       (1 + epoch * t | user_id)),
  data = bysubj_rev_df,
  family = student(),
  backend = "cmdstanr",
  chains = 4,
  iter = 2000,
  prior = c(
    set_prior("normal(0, 1)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 2)", class = "sd", group = "user_id"), # Weakly informative prior for random effects sd
    set_prior("lkj(1)", class = "cor", group = "user_id") # Weakly informative prior for random effects correlations
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  threads = threading(2)
)

pp_check(fit_rev) + xlim(-2, 2)

summary(fit_rev)
marginal_effects(fit_rev, "t")
marginal_effects(fit_rev, "epoch")
marginal_effects(fit_rev, "t:epoch")

performance::r2_bayes(fit_rev)

(loo_fit_rev <- loo(fit_rev))
# This model has no divergent transitions and only 2 'bad' k diagnostic values.

fit_rev_no_t <- update(
  fit_rev,
  formula = bf(happiness ~ 1 + epoch + (1 + epoch | user_id)) # remove the t term
)

loo_fit_rev_no_t <- loo(fit_rev_no_t)

loo_compare(loo_fit_rev, loo_fit_rev_no_t)
#             elpd_diff se_diff
# fit_rev         0.0       0.0 
# fit_rev_no_t -295.8      53.4 

report::report_effectsize(fit_rev, effectsize_method = "basic")
# very small (Std. beta = 0.00, 95% CI [0.00, 0.00])
# very small (Std. beta = -0.07, 95% CI [-0.09, -0.05])
# very small (Std. beta = 0.08, 95% CI [0.05, 0.11])
# very small (Std. beta = -0.01, 95% CI [-0.03, 7.26e-03])
# very small (Std. beta = 0.11, 95% CI [0.09, 0.12])

# Random slopes

# Assuming that the naming convention for brms is such that the random slopes for t are named "r_user_idt" or similar
random_slope_t_pattern <- "^r_user_id\\[.*,t\\]$"

# Extract only the random slopes for `t` for each user
random_slopes_t <- random_effects_samples[ , grepl(random_slope_t_pattern, colnames(random_effects_samples))]

# Now we should have one column for each user's random effect for t.
# The number of columns should match the number of users
if (ncol(random_slopes_t) != length(unique(bysubj_rev_df$user_id))) {
  stop("Unexpected number of random slopes for t extracted")
}

# Get the posterior mean of the fixed effect for `t`
fixed_effect_t <- apply(fixed_effects_samples, 2, mean)["b_t"]

# Get the mean random slope for `t` across all posterior samples for each user
mean_random_slopes_t <- apply(random_slopes_t, 2, mean)

# Total slope for `t` for each user is the sum of the fixed effect and their random effect
mean_total_slopes_t <- fixed_effect_t + mean_random_slopes_t
mean_total_slopes_t <- as.vector(mean_total_slopes_t)

# Check the length of the final total slopes vector
length(mean_total_slopes_t)  # Should now be 215

hist(mean_total_slopes_t)

plot(total_slopes_t, mean_total_slopes_t)
cor(total_slopes_t, mean_total_slopes_t)

mean(mean_total_slopes_t > 0)
# [1] 0.7581395
# The proportion of positive slopes from brms is the same as that found 
# with lmer.

# ------------------

# Fit the second model with the 'brm' function from the 'brms' package
# The model is parametrized so as to test directly a difference in slope and
# a difference in intercept between the two regression lines before and after
# the break point which distinguishes between the two epochs.
# IMPORTANT: This model can only be used for specific tests on intercept 
# difference and slope differences between the two segments. It is better to use 
# the previous model for the general test of the effect of trial on momentary 
# happiness.
fit2_rev <- brm(
  bf(
    happiness ~ 0 + Intercept + epoch:Intercept + t + epoch:t + 
      (0 + t + epoch:t | user_id)
    # The intercept is the base intercept for t < 0.
    # The epoch:intercept is the difference in intercept for t >= 0.
    # The t is the slope for t < 0.
    # The epoch:t is the difference in slope for t >= 0.
  ),
  data = bysubj_rev_df,
  family = student(),
  backend = "cmdstanr",
  cores = 8, 
  chains = 4,
  iter = 5000,
  prior = c(
    set_prior("normal(0, 1)", class = "b"), # Weakly informative prior for fixed effects
    set_prior("cauchy(0, 1)", class = "sd"), # Weakly informative prior for random effects sd
    set_prior("lkj(1)", class = "cor") # Weakly informative prior for random effects correlations
  ),
  threads = threading(2),
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

pp_check(fit2_rev) + xlim(-2, 2)

traceplot(
  fit2_rev, 
  pars = c("muc", "mur", "muu", "Omega", "sigma", "lp__")
)

summary(fit2_rev)
conditional_effects(fit2_rev, "t:epoch")
conditional_effects(fit2_rev, "t")

performance::r2_bayes(fit2_rev)



# eof ----


