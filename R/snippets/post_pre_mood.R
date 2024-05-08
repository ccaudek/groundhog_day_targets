#' Examine the mood difference post-mood - pre_mood as a function 
#' of is_reversal.


tar_load(params_happiness_clean_df)

# Calcola la differenza
params_happiness_clean_df <- params_happiness_clean_df %>%
  mutate(mood_diff = mood_post - mood_pre)

# Standardizzazione per user_id
params_happiness_clean_df <- params_happiness_clean_df |> 
  group_by(user_id) |> 
  mutate(
    mean_diff = mean(mood_diff, na.rm = TRUE),
    sd_diff = sd(mood_diff, na.rm = TRUE),
    standardized_mood_diff = (mood_diff - mean_diff) / sd_diff,
    mean_mood_pre = mean(mood_pre, na.rm = TRUE),
    sd_mood_pre = sd(mood_pre, na.rm = TRUE),
    standardized_mood_pre = (mood_pre - mean_mood_pre) / sd_mood_pre
  ) |> 
  ungroup()


mydat <- params_happiness_clean_df |> 
  group_by(user_id, ema_number, is_reversal) |> 
  summarize(
    zdelta_mood = mean(standardized_mood_diff, na.rm = TRUE),
    moodc_pre = mean(standardized_mood_pre, na.rm = TRUE)
  ) |> 
  ungroup()

mydat$is_reversal <- factor(mydat$is_reversal)
mydat$zt <- as.vector(scale(mydat$ema_number))


get_prior(
  zdelta_mood ~ 1 + zt + moodc_pre * is_reversal + (1 + zt | user_id),
  data = mydat
)

my_priors <- c(
  set_prior("normal(0, 1)", class = "b", coef = "moodc_pre"),
  set_prior("normal(0, 1)", class = "b", coef = "zt"),
  set_prior("normal(0, 1)", class = "b", coef = "is_reversal1")
)

mod <- brm(
  zdelta_mood ~ 1 + zt + moodc_pre * is_reversal + (1 + zt | user_id),
  prior = my_priors, 
  backend = "cmdstanr",
  family = student(),
  data = mydat,
  # control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
  chains = 4, cores = 8, threads = threading(2),
  silent = 2
)
pp_check(mod)

loo_mod <- loo(mod)
plot(loo_mod)

summary(mod)

performance::r2_bayes(mod)
# Conditional R2: 0.696 (95% CI [0.688, 0.704])
#    Marginal R2: 0.692 (95% CI [0.684, 0.699])

conditional_effects(mod, "moodc_pre")
conditional_effects(mod, "is_reversal")

effectsize::standardize_parameters(
  mod
)
# Standardization method: refit

# Component   |                Parameter | Std. Median |         95% CI
# ---------------------------------------------------------------------
# conditional |              b_Intercept |        0.06 | [ 0.02,  0.10]
# conditional |                     b_zt |        0.02 | [ 0.00,  0.04]
# conditional |              b_moodc_pre |       -0.84 | [-0.88, -0.80]
# conditional |           b_is_reversal1 |       -0.10 | [-0.14, -0.05]
# conditional | b_moodc_pre:is_reversal1 |       -0.09 | [-0.13, -0.04]
# sigma       |                    sigma |        0.34 | [ 0.31,  0.37]

my_priors_bis <- c(
  set_prior("normal(0, 1)", class = "b", coef = "moodc_pre"),
  set_prior("normal(0, 1)", class = "b", coef = "zt")
)

mod_bis <- brm(
  zdelta_mood ~ 1 + zt + moodc_pre + (1 + zt | user_id),
  prior = my_priors_bis, 
  backend = "cmdstanr",
  family = student(),
  data = mydat,
  # control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
  chains = 4, cores = 8, threads = threading(2),
  silent = 2
)


loo_mod_bis <- loo(mod_bis)
plot(loo_mod_bis)

loo_compare(loo_mod_bis, loo_mod)        
#       elpd_diff se_diff
# mod       0.0       0.0  
# mod_bis -14.1       6.1  


foo <- emmeans(mod1, ~ time + is_reversal)




