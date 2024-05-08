# Groundhog Day project
# PUUI: 9811beb3-2dc5-481f-8b0f-2e8c12936507

# This script is used to test the targets pipeline.

# Libraries

suppressPackageStartupMessages({
  library("here")
  library("tidyverse")
  library("mice")
  library("lme4")
  library("brms")
  library("bayesplot")
  library("effectsize")
  library("scales")
  library("sjstats")
  library("sjPlot")
  library("sjmisc")
})


source(here::here("R", "funs", "funs_instant_mood.R"))
source(here::here("R", "funs", "funs_quest.R"))
# Meta-functions
source(here::here("R", "functions.R"))


# Get data ---------------------------------------------------------------------

# Get raw data and remove bad subjects.
file_path <- here::here("data", "prep", "groundhog_all_clean.RDS")
MAX_EMA_NUMBER <- 14 # Do not use EMA sessions greater than 14.

df <- get_data(file_path, MAX_EMA_NUMBER)


# Get parameters of the momentary happiness model ------------------------------

# Get list of unique user_ids
user_id_codes <- unique(df$user_id)
# Get data frame with the parameters for each ema session and user
params_happiness_df <- get_params_happiness_model(df, user_id_codes)


# Clean parameters -------------------------------------------------------------

params_happiness_clean_df <- clean_params_happiness_model(params_happiness_df)


# alpha analysis ---------------------------------------------------------------

params_happiness_clean_df$alpha |> hist()

mod_alpha <- brm(
  alpha ~ environment +
    (environment | user_id / ema_number),
  family = asym_laplace(),
  algorithm = "meanfield",
  data = params_happiness_clean_df,
  file = here::here("R", "brms", "mod_alpha.RDS"),
  refresh = 0
)

pp_check(mod_alpha) + xlim(0, 1)
summary(mod_alpha)
bayes_R2(mod_alpha)
conditional_effects(mod_alpha, "environment")


# Mood Difference --------------------------------------------------------------

tar_load(params_happiness_clean_df)

params_happiness_clean_df <- params_happiness_clean_df |> 
  group_by(user_id, environment, ema_number) |> 
  mutate(
    mood_dif = mood_post - mood_pre
  )

params_happiness_clean_df |> 
  group_by(user_id, environment) |> 
  summarize(
    mood_pre = mean(mood_pre, trim = 0.1),
    mood_post = mean(mood_post, trim = 0.1),
    mood_dif = mean(mood_dif, trim = 0.1)
  ) |> 
  group_by(environment) |> 
  summarize(
    mood_pre = mean(mood_pre, trim = 0.1),
    mood_post = mean(mood_post, trim = 0.1),
    mood_dif = mean(mood_dif, trim = 0.1)
  )

params_happiness_clean_df <- params_happiness_clean_df |> 
  group_by(user_id) |> 
  mutate(
    mood_pre_cw = mood_pre - mean(mood_pre, trim = 0.1)
  )
  
mod3 <- brm(
  mood_dif ~ mood_pre_cw + environment +
    (mood_pre_cw + environment | user_id / ema_number),
  family = student(),
  algorithm = "meanfield",
  data = params_happiness_clean_df,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

pp_check(mod3) + xlim(-100, 100)
summary(mod3)
conditional_effects(mod3, "ema_number:environment")
bayes_R2(mod3)


mod_mood_2 <- brm(
  mood_dif ~ mood_pre_cw + ema_number +
    environment * (w1 + w2 + w3 + w4 + w5 + gamma) +
    (mood_pre_cw + ema_number + environment  | user_id ),
  family = student(),
  algorithm = "meanfield",
  cores = 4,
  data = params_happiness_clean_df,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  refresh = 0
)


tar_load(params_happiness_df)

# Unnest the mle_params column
expanded_df <- params_happiness_df %>%
  unnest_wider(mle_params, names_sep = "_")

expanded_df <- expanded_df %>% 
  rename(
    w0 = mle_params_1,
    w_outcome = mle_params_2,
    w_stimulus = mle_params_3,
    w_rpe = mle_params_4,
    w_moodpre = mle_params_5,
    w_control = mle_params_6,
    w_trial = mle_params_7,
    w_gamma = mle_params_8
  )

glimpse(expanded_df)


outlier_cols1 <- c(
  "w0", "w_outcome", "w_stimulus", "w_rpe", "w_moodpre", "w_control", 
  "w_trial", "w_gamma"
  )

params_happiness_df <- 
  impute_outliers_params_happiness(expanded_df, outlier_cols1)

# Step 2: Detect and impute outliers for mood and control variables
outlier_cols2 <- c("mood_pre", "control", "mood_post")
params_happiness_df <- 
  impute_outliers_mood(params_happiness_df, outlier_cols2)

# Step 3: Additional transformations
params_happiness_df <- params_happiness_df %>%
  group_by(user_id) %>%
  mutate(
    mood_dif = mood_post - mood_pre,
    mood_pre_cw = mood_pre - mean(mood_pre, trim = 0.1),
    environment = ifelse(is_reversal == 1, "Volatile", "Stable")
  ) %>%
  ungroup()






