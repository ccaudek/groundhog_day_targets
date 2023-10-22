# _targets.R #

# Created by use_targets()
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   <https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline> # nolint

# Load packages required to define the pipeline:
library("here")
library("targets")
library("tarchetypes")

library("dplyr")
library("mice")
library("brms")
library("bayesplot")
library("effectsize")
library("scales")
library("sjstats")
library("sjPlot")
library("sjmisc")


# Set target options:
tar_option_set(
  # packages that your targets need to run
  packages = c(
    "here", "tidyverse", "mice", "brms", "quarto", "bayesplot", 
    "effectsize", "sjstats", "sjPlot", "sjmisc", "emmeans"), 
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
source(here::here("R", "funs", "funs_instant_mood.R"))
source(here::here("R", "funs", "funs_quest.R"))
source(here::here("R", "functions.R"))
# source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  
  # Get raw PRL data.
  tar_target(
    raw_df, 
    here::here("data", "prep", "groundhog_all_clean.RDS"), format = "file"
  ),
  
  # Clean PRL data. raw_df is the input to get_data() and prl_df is the output.
  tar_target(prl_df, get_data(raw_df)),
  
  # Get parameters of the momentary happiness model.
  tar_target(
    params_happiness_df, 
    get_params_happiness_model(prl_df, unique(prl_df$user_id))
  ),

  # Clean parameters of the momentary happiness model.
  tar_target(
    params_happiness_clean_df, 
    clean_params_happiness_model(params_happiness_df)
  ),
  
  # brms model for alpha and volatility
  tar_target(
    brms_fitted_mod_alpha, 
    brms_mod_alpha(params_happiness_clean_df)
  ),
  
  # brms model for mood_dif 
  tar_target(
    brms_fitted_mod_mood_1, 
    brms_mod_mood_1(params_happiness_clean_df)
  ),
  
  # brms model for mood_dif with also happiness parameters
  # tar_target(
  #   brms_fitted_mod_mood_2, 
  #   brms_mod_mood_2(params_happiness_clean_df)
  # ),
  
  # Create report
  tar_quarto(
    name = report,
    path = here::here("doc/groundhog_day_report.qmd")
  )
  
# close list  
)





