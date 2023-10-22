
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
