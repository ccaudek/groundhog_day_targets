
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
  
  # The data have already been cleaned so that each participant completed
  # at least 4 EMA sessions.
  
  # file_path <- here::here("data", "prep", "groundhog_all_clean.RDS")
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
  
  # Compute ema_number and standardize mood and control by user_id.
  d4 <- d3 |> 
    group_by(user_id) |> 
    arrange(date) |>  
    mutate(
      ema_number = dense_rank(date),
      zmoodpre = as.vector(scale(mood_pre, center = TRUE, scale = TRUE)),
      zcontrol = as.vector(scale(control, center = TRUE, scale = TRUE)),
      zim = as.vector(scale(instant_mood, center = TRUE, scale = TRUE))
    ) |> 
    ungroup()
  
  # table(d4$ema_number)
  #    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
  # 6450 6450 6450 6450 6090 5820 5550 5250 4890 4170 3600 2760 1890  390   60   30   30 
  
  d5 <- d4 |> 
    dplyr::filter(ema_number < 13)
  
  d6 <- d5 |> 
    dplyr::select(
      instant_mood, feedback, trial, user_id, mood_pre, control, mood_post,
      is_reversal, is_target_chosen, ema_number, zim, zmoodpre, zcontrol, 
      date
    )

  d6
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
get_params_happiness_model <- function(dz_clean, user_id_codes) {
  
  # Apply the function process_user() to each user_id.
  results_list <- lapply(user_id_codes, process_user, dz_clean)
  
  # Bind all data frames together into a single data frame.
  all_results_df <- bind_rows(results_list)
  
  # Add mood_pre, mood_post, control.
  bysubj_mood_df <- dz_clean %>%
    group_by(user_id, ema_number) %>%
    summarize(
      mood_pre = mean(mood_pre),
      mood_post = mean(mood_post),
      control = mean(control),
      zmoodpre = mean(zmoodpre),
      zcontrol = mean(zcontrol)
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
  results_df$control <- remove_outliers(results_df$control)
  
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
    iter = 5000,
    warmup = 1000,
    cores = 8,
    data = params_happiness_clean_df,
    control = list(adapt_delta = 0.99, max_treedepth = 20),
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


# brms_mod_mood_3() ------------------------------------------------------------

brms_mod_mood_3 <- function(params_happiness_clean_df) {
  # Create the quadratic term
  params_happiness_clean_df$ema_number_sq <-
    params_happiness_clean_df$ema_number^2
  params_happiness_clean_df$ema_number_cu <-
    params_happiness_clean_df$ema_number^3

  # Center and scale the linear and quadratic terms
  params_happiness_clean_df$ema_number <-
    scale(params_happiness_clean_df$ema_number)
  params_happiness_clean_df$ema_number_sq <-
    scale(params_happiness_clean_df$ema_number_sq)
  params_happiness_clean_df$ema_number_cu <-
    scale(params_happiness_clean_df$ema_number_cu)

  # Fit the model
  mod_mood_3 <- brm(
    mood_dif ~ mood_pre_cw +
      environment * (ema_number + ema_number_sq + ema_number_cu) +
      (mood_pre_cw +
        environment * (ema_number + ema_number_sq + ema_number_cu) | user_id
      ),
    data = params_happiness_clean_df,
    family = student(),
    backend = "cmdstanr",
    cores = 4,
    refresh = 0
  )
}


# get_plot_mood_dif_ema_number() -----------------------------------------------

get_plot_mood_dif_ema_number <- function(
    params_happiness_clean_df, brms_fitted_mod_mood_3) {
  
  ## For debugging purposes
  # params_happiness_clean_df |>
  #   group_by(environment, ema_number) |>
  #   summarize(
  #     md = mean(mood_dif),
  #     se = sqrt(var(mood_dif) / n()),
  #     .groups = 'drop'  # This argument drops the grouping structure afterwards
  #   ) |>
  #   ggplot(aes(x = ema_number, y = md, color = environment, group = environment)) +
  #   geom_line() +
  #   geom_ribbon(aes(ymin = md - se, ymax = md + se, fill = environment), alpha = 0.2) +
  #   facet_wrap(~ environment)

  # Create the quadratic term
  params_happiness_clean_df$ema_number_sq <-
    params_happiness_clean_df$ema_number^2
  params_happiness_clean_df$ema_number_cu <-
    params_happiness_clean_df$ema_number^3

  # Center and scale the linear and quadratic terms
  params_happiness_clean_df$ema_number <-
    scale(params_happiness_clean_df$ema_number)
  params_happiness_clean_df$ema_number_sq <-
    scale(params_happiness_clean_df$ema_number_sq)
  params_happiness_clean_df$ema_number_cu <-
    scale(params_happiness_clean_df$ema_number_cu)

  # Create a new data frame for predictions for each level of 'environment'
  env_levels <- unique(params_happiness_clean_df$environment)
  prediction_list <- lapply(env_levels, function(env) {
    new_data <- data.frame(
      ema_number = scale(seq(min(params_happiness_clean_df$ema_number), max(params_happiness_clean_df$ema_number), length.out = 100)),
      mood_pre_cw = mean(params_happiness_clean_df$mood_pre_cw, na.rm = TRUE),
      environment = env,
      user_id = unique(params_happiness_clean_df$user_id)[1] # Replace with an actual user_id if needed
    )

    # Add the scaled quadratic and cubic terms
    new_data$ema_number_sq <- scale(new_data$ema_number^2)
    new_data$ema_number_cu <- scale(new_data$ema_number^3)

    # Generate predictions
    preds <- predict(brms_fitted_mod_mood_3, newdata = new_data, summary = FALSE, re.form = NA)

    # Summarize predictions
    summary_preds <- apply(preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    summary_preds_df <- as.data.frame(t(summary_preds))
    colnames(summary_preds_df) <- c("Q2.5", "Estimate", "Q97.5")

    # Combine with new_data and return
    cbind(new_data, setNames(summary_preds_df, paste0("pred_", colnames(summary_preds_df))))
  })

  # Combine all predictions into a single data frame
  all_predictions <- do.call(rbind, prediction_list)

  # Okabe-Ito colors
  okabe_ito_colors <- c("Stable" = "#56B4E9", "Volatile" = "#D55E00")

  # Plot
  plot_mood_dif_ema_number <- ggplot(all_predictions, aes(x = ema_number)) +
    geom_line(aes(y = pred_Estimate, color = environment)) +
    geom_ribbon(
      aes(ymin = pred_Q2.5, ymax = pred_Q97.5, fill = environment),
      alpha = 0.4
    ) +
    scale_color_manual(values = okabe_ito_colors) +
    scale_fill_manual(values = okabe_ito_colors) +
    labs(
      title = "Combined Effect of Linear, Quadratic, and Cubic \nComponents by Environment",
      x = "EMA Number (scaled)",
      y = "Predicted Mood Difference"
    )

  plot_mood_dif_ema_number
}

