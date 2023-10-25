
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


