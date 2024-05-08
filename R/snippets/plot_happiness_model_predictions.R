# Plot the average value of momentary mood as a function of trial, for 
# reversal and no-reversal sessions. Adds the predictions of the momentary
# happiness model.

tar_load(params_happiness_clean_df)

# Plot happiness and model's predictions

# Generate a plot of momentary happiness as a function of trial and 
# is_reversal. Add to the plot the predictions of the momentary happiness
# model.

# Calculate mean and standard error of the mean (SEM)
for_plot_df <- params_happiness_clean_df |> 
  group_by(is_reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(predicted_happiness, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(predicted_happiness, na.rm = TRUE) / sqrt(n())
  ) |> 
  ungroup()
for_plot_df$is_reversal <- as.factor(for_plot_df$is_reversal)

# Create the plot
for_plot_df |> 
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = viridis(3)[1], linewidth = 2) +  # Line for h
  geom_ribbon(aes(ymin = h - h_sem, ymax = h + h_sem), alpha = 0.2) +  # Ribbon for h
  geom_line(aes(y=h_hat), color = viridis(3)[2], linewidth = 1.0) +  # Line for h_hat
  #geom_ribbon(aes(ymin = h_hat - h_hat_sem, ymax = h_hat + h_hat_sem), alpha = 0.2, fill = "red") +  # Ribbon for h_hat
  facet_wrap(~ is_reversal) +
  theme_default()  # Apply the bayesplot theme


with(
  params_happiness_clean_df,
  cor(happiness, predicted_happiness
)








  
  
)