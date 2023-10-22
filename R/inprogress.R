
# Counting the number of unique dates for each user_id
d1 <- readRDS("data/prep/groundhog_all_clean.RDS")


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

unique_dates_per_user <- d3 %>% 
  group_by(user_id) %>% 
  summarise(num_unique_dates = n_distinct(date))

table(unique_dates_per_user$num_unique_dates)


prl_df_ranked <- d3 %>%
  group_by(user_id) %>% 
  arrange(date) %>% 
  mutate(date_rank = dense_rank(date)) %>% 
  ungroup()

table(prl_df_ranked$date_rank)

params_happiness_clean_df |> 
  group_by(environment, ema_number) |> 
  summarize(
    md = mean(mood_dif),
    se = sqrt(var(mood_dif) / n())
  ) |> 
  ggplot(aes(x=ema_number, y = md, color = environment)) +
  geom_line() 
  #facet_wrap(~ environment)

# Create the quadratic term
params_happiness_clean_df$ema_number_sq <- params_happiness_clean_df$ema_number^2
params_happiness_clean_df$ema_number_cu <- params_happiness_clean_df$ema_number^3

# Center and scale the linear and quadratic terms
params_happiness_clean_df$ema_number <- scale(params_happiness_clean_df$ema_number)
params_happiness_clean_df$ema_number_sq <- scale(params_happiness_clean_df$ema_number_sq)
params_happiness_clean_df$ema_number_cu <- scale(params_happiness_clean_df$ema_number_cu)


mod_mood_2 <- brm(
  mood_dif ~ mood_pre_cw + environment + 
    ema_number + ema_number_sq + ema_number_cu +
    (mood_pre_cw + environment + 
       ema_number + ema_number_sq + ema_number_cu | user_id),
  family = student(),
  algorithm = "meanfield",
  # iter = 5000,
  # warmup = 1000,
  cores = 8,
  data = params_happiness_clean_df
  # control = list(adapt_delta = 0.99, max_treedepth = 20),
)
pp_check(mod_mood_2)
summary(mod_mood_2)

conditional_effects(mod_mood_2, "ema_number")
conditional_effects(mod_mood_2, "ema_number_sq")
conditional_effects(mod_mood_2, "ema_number_cu")





# Fit the model
mod_mood_2 <- brm(
  mood_dif ~ mood_pre_cw + environment * (ema_number + ema_number_sq + ema_number_cu) +
    (mood_pre_cw + environment * (ema_number + ema_number_sq + ema_number_cu) | user_id),
  data = params_happiness_clean_df,
  family = student(),
  algorithm = "meanfield",
  cores = 8
)

# Create a new data frame for predictions for each level of 'environment'
env_levels <- unique(params_happiness_clean_df$environment)
prediction_list <- lapply(env_levels, function(env) {
  new_data <- data.frame(
    ema_number = scale(seq(min(params_happiness_clean_df$ema_number), max(params_happiness_clean_df$ema_number), length.out = 100)),
    mood_pre_cw = mean(params_happiness_clean_df$mood_pre_cw, na.rm = TRUE),
    environment = env,
    user_id = unique(params_happiness_clean_df$user_id)[1]  # Replace with an actual user_id if needed
  )
  
  # Add the scaled quadratic and cubic terms
  new_data$ema_number_sq <- scale(new_data$ema_number^2)
  new_data$ema_number_cu <- scale(new_data$ema_number^3)
  
  # Generate predictions
  preds <- predict(mod_mood_2, newdata = new_data, summary = FALSE, re.form = NA)
  
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
ggplot(all_predictions, aes(x = ema_number)) +
  geom_line(aes(y = pred_Estimate, color = environment)) +
  geom_ribbon(aes(ymin = pred_Q2.5, ymax = pred_Q97.5, fill = environment), alpha = 0.4) +
  scale_color_manual(values = okabe_ito_colors) +
  scale_fill_manual(values = okabe_ito_colors) +
  labs(
    title = "Combined Effect of Linear, Quadratic, and Cubic \nComponents by Environment",
    x = "EMA Number (scaled)",
    y = "Predicted Mood Difference"
  )

