# _targets.R #

# Created by use_targets()
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   <https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline> # nolint

# Load packages required to define the pipeline:
suppressPackageStartupMessages({
  library("here")
  library("targets")
  library("tarchetypes")
  library("tidyverse")
  library("brms")
  library("mice")
  library("parallel")
  library("viridis")
  library("bayesplot") 
})


# Set target options:
tar_option_set(
  # packages that your targets need to run
  packages = c(
    "here", "tarchetypes", "tidyverse", "mice", "brms", "quarto", "bayesplot", 
    "effectsize", "sjstats", "sjPlot", "sjmisc", "emmeans", "parallel", 
    "viridis", "bayesplot"), 
  format = "rds", # default storage format
  seed = 12345
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source(here::here("R", "funs", "funs_mood_modeling.R"))
tar_source(here::here("R", "funs", "funs_brms.R"))
tar_source(here::here("R", "funs", "funs_quest.R"))


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
  # nrow(params_happiness_df) = 2131
  # names(params_happiness_df)
  # [1] "w0"          "w1"          "w2"          "w3"          "w4"         
  # [6] "w5"          "gamma"       "is_reversal" "ema_number"  "user_id"    
  # [11] "alpha"       "mood_pre"    "mood_post"   "control"     "zmoodpre"   
  # [16] "zcontrol"  
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
  # tar_target(
  #   brms_fitted_mod_alpha,
  #   brms_mod_alpha(params_happiness_clean_df)
  # ),

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




