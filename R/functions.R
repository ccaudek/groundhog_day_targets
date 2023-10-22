
remove_outliers <- function(my_vector) {
  # Calculate mean and standard deviation
  mean_val <- mean(my_vector, na.rm = TRUE)
  std_dev <- sd(my_vector, na.rm = TRUE)
  # Calculate Z-scores
  z_scores <- (my_vector - mean_val) / std_dev
  # Replace outliers with NA
  my_vector[abs(z_scores) > 3.5] <- NA
  
  my_vector
}


# Import PRL data --------------------------------------------------------------
#'
#' @description
#' Read pre-processed complete raw data and remove bad subjects.
#' @param file_path The path to the pre-processed RDS file with all raw data.
#' @returns A dataframe.
#' 
get_data <- function(file_path) {
  
  d1 <- readRDS(file_path) 
  
  d1$user_id <- as.numeric(d1$user_id)
  
  # Remove user_id == NA.
  d2 <- d1 |>
    dplyr::filter(!is.na(user_id))
  
  # user_id with convergence problems.
  bad_ids <- c(
    3338029881, 3665345709, 3248648540, 3668615349, 3456996255, 
    3475648095, 3497824506, 3888185906, 3667000898
  )
  
  # Remove bad user_id.
  d3 <- d2[!(d2$user_id %in% bad_ids), ]
  
  # The data have already been cleaned so that each participant completed
  # at least 4 EMA sessions.
  
  d4 <- d3 |> 
    group_by(user_id) |> 
    arrange(date) |>  
    mutate(ema_number = dense_rank(date)) |> 
    ungroup()
  
  # table(d4$ema_number)
  #    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
  # 6450 6450 6450 6450 6090 5820 5550 5250 4890 4170 3600 2760 1890  390   60   30   30 
  
  d5 <- d4 |> 
    dplyr::filter(ema_number < 13)

  d5
}


# get_params_happiness_model() -------------------------------------------------
#'
#' @description
#' Get a data frame with the parameters of the subjective momentary
#' happiness model for each participant and each ema session, together 
#' with control, mood_pre, and mood_post.
#' 
#' @param d_clean A data frame.
#' @returns A data frame.
#' 
get_params_happiness_model <- function(d_clean, user_id_codes) {
  
  # Standardize momentary mood by user_id.
  dz_clean <- d_clean %>%
    group_by(user_id) %>%
    mutate(
      zim = as.vector(scale(instant_mood, center = TRUE, scale = TRUE))
    ) %>%
    ungroup()
  
  # Apply the function process_user() to each user_id.
  results_list <- lapply(user_id_codes, process_user, dz_clean = dz_clean)
  
  # Bind all data frames together into a single data frame.
  all_results_df <- bind_rows(results_list)
  
  # Add mood_pre, mood_post, control.
  bysubj_mood_df <- d_clean %>%
    group_by(user_id, ema_number) %>%
    summarize(
      mood_pre = mean(mood_pre),
      mood_post = mean(mood_post),
      control = mean(control)
    ) %>%
    ungroup()
  
  # Final dataframe.
  results_df <- left_join(
    all_results_df, bysubj_mood_df,
    by = c("user_id", "ema_number")
  )
  
  return(results_df)
}


# clean_params_happiness_model() -----------------------------------------------

clean_params_happiness_model <- function(results_df) {
  
  results_df$w1 <- remove_outliers(results_df$w1)
  results_df$w2 <- remove_outliers(results_df$w2)
  results_df$w3 <- remove_outliers(results_df$w3)
  results_df$w4 <- remove_outliers(results_df$w4)
  results_df$w5 <- remove_outliers(results_df$w5)
  results_df$gamma <- remove_outliers(results_df$gamma)
  results_df$mood_pre <- remove_outliers(results_df$mood_pre)
  results_df$mood_post <- remove_outliers(results_df$mood_post)
  
  # Replace with NA the values alpha = 0 ad alpha = 1.
  results_df$alpha <- ifelse(
    results_df$alpha < 0.00001 | results_df$alpha > 0.99999, 
    NA, results_df$alpha
  )
  
  # Multiple imputation
  imputed_cart <- complete(mice(results_df, method = "cart"))
  
  # On the cleaned and imputed data, compute the difference 
  # mood_post - mood_pre.
  bysubj_df <- imputed_cart |>
    group_by(user_id) |> 
    mutate(
      mood_dif = mood_post - mood_pre,
      mood_pre_cw = mood_pre - mean(mood_pre, trim = 0.1)
    ) 
  
  bysubj_df$environment <- ifelse(
    bysubj_df$is_reversal == 1, "Volatile", "Stable"
  ) 
  
  bysubj_df
}


# brms_mod_alpha() -------------------------------------------------------------

brms_mod_alpha <- function(params_happiness_clean_df) {
  
  mod_alpha <- brm(
    alpha ~ environment +
      (environment | user_id / ema_number),
    family = asym_laplace(),
    algorithm = "meanfield",
    data = params_happiness_clean_df,
    refresh = 0
  )
}


# brms_mod_mood_1() ------------------------------------------------------------

brms_mod_mood_1 <- function(params_happiness_clean_df) {

  mod_mood_1 <- brm(
    mood_dif ~ mood_pre_cw + environment +
      (mood_pre_cw + environment | user_id / ema_number),
    family = student(),
    backend = "cmdstanr",
    # iter = 5000,
    # warmup = 1000,
    cores = 8,
    data = params_happiness_clean_df,
    # control = list(adapt_delta = 0.99, max_treedepth = 20),
    refresh = 0
  )
}


# brms_mod_mood_2() ------------------------------------------------------------

brms_mod_mood_2 <- function(params_happiness_clean_df) {

  mod_mood_2 <- brm(
    mood_dif ~ control + mood_pre_cw + 
      environment * (w1 + w2 + w3 + w4 + w5 + gamma) +
      (control + mood_pre_cw + environment | user_id / ema_number),
    family = student(),
    backend = "cmdstanr",
    # iter = 5000,
    # warmup = 1000,
    cores = 8,
    data = params_happiness_clean_df,
    # control = list(adapt_delta = 0.99, max_treedepth = 20),
    refresh = 0
  )
}



