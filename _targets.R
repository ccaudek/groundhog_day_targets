# _targets.R #

# Created by use_targets()
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   <https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline> # nolint

# Loading package namespaces:
suppressPackageStartupMessages({
  library("here")
  library("targets")
  library("tarchetypes")
  library("tidyverse")
  library("tidyr")
  library("purrr")
  library("cmdstanr")
  library("brms")
  library("posterior")
  library("loo")
  library("lme4")
  library("mice")
  library("parallel")
  library("emmeans")
  library("quarto")
  library("bayesplot") 
  library("gridExtra")
  library("ggdist")
  library("viridis")
})

options(pillar.neg=FALSE)
theme_set(bayesplot::theme_default(base_family = "sans"))

# Set random seed for reproducibility:
SEED <- 48927 

# Set target options:
tar_option_set(
  # packages that your targets need to run
  packages = c(
    "here", "targets", "tarchetypes", "tidyverse", "tidyr", "purrr",
    "cmdstanr", "brms", "posterior", "loo", "mice", "parallel", "emmeans",
    "quarto", "bayesplot", "gridExtra", "ggdist", "effectsize",
    "sjstats", "sjPlot", "sjmisc", "viridis"
  ),
  format = "rds", # default storage format
  seed = SEED
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Sourcing helper files:
list.files(
  here::here("R", "funs"),
  pattern = "\\.R$", recursive = TRUE,
  full.names = TRUE
) %>%
  walk(source)


# Replace the target list below with your own:
list(
  
  # Get raw PRL data.
  tar_target(
    raw_df,
    here::here("data", "prep", "groundhog_all_clean.RDS"), format = "file"
  ),

  # Clean PRL data. raw_df is the input to get_data() and prl_df is the output.
  # prl_df has 63930 rows.
  tar_target(prl_df, get_data(raw_df)),

  # Estimate the parameters of the momentary happiness model.
  tar_target(
    params_happiness_df,
    get_params_happiness_model(prl_df, unique(prl_df$user_id))
  ),

  # Clean parameters of the momentary happiness model.
  tar_target(
    params_happiness_clean_df,
    clean_params_happiness_model(params_happiness_df)
  )

  # brms model for alpha and volatility
  # the distribution of alpha is bimodal. No congergence!
  # tar_target(
  #   brms_fitted_mod_alpha,
  #   brms_mod_alpha(params_happiness_clean_df)
  # )

  # brms model for mood_dif
  # tar_target(
  #   brms_fitted_mod_mood_1,
  #   brms_mod_mood_1(params_happiness_clean_df)
  # ),
  
  # brms model for mood_dif with also happiness parameters
  # tar_target(
  #   brms_fitted_mod_mood_2, 
  #   brms_mod_mood_2(params_happiness_clean_df)
  # ),
  
  # Fit the brms model with mood_dif as a function of time x environment.
  # mood_pre is used as a covariate, and time is coded in terms of the linear, 
  # quadratic, and cubic components of ema_number.
  # tar_target(
  #   brms_fitted_mod_mood_3,
  #   brms_mod_mood_3(params_happiness_clean_df)
  # ),
  
  # Compute predicted values of model brms_fitted_mod_mood_3 and 
  # generate figure.
  # tar_target(
  #   plot_mood_dif_ema_number,
  #   get_plot_mood_dif_ema_number(
  #     params_happiness_clean_df, brms_fitted_mod_mood_3
  #   )
  # ),
  # 
  # Create report.
  # tar_quarto(
  #   name = report,
  #   path = here::here("doc", "groundhog_day_report.qmd")
  # )
  
# close list  
)



