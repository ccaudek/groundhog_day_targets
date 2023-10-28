

# Calculate mean and standard error of the mean (SEM)
foo <- final_df |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1, na.rm = TRUE),
    h_hat = mean(happiness_hat, trim = 0.1, na.rm = TRUE),
    h_sem = sd(happiness, na.rm = TRUE) / sqrt(n()),
    h_hat_sem = sd(happiness_hat, na.rm = TRUE) / sqrt(n())
  ) |> 
  ungroup()
foo$reversal <- as.factor(foo$reversal)

# Create the plot
foo |> 
  ggplot(aes(x=trial)) +
  geom_line(aes(y=h), color = "blue") +  # Line for h
  geom_ribbon(aes(ymin = h - h_sem, ymax = h + h_sem), alpha = 0.2) +  # Ribbon for h
  geom_line(aes(y=h_hat), color = "red") +  # Line for h_hat
  geom_ribbon(aes(ymin = h_hat - h_hat_sem, ymax = h_hat + h_hat_sem), alpha = 0.2, fill = "red") +  # Ribbon for h_hat
  facet_wrap( ~ reversal)



glimpse(params_happiness_clean_df)

params_happiness_clean_df$session <- 
  as.vector(params_happiness_clean_df$ema_number)

m1 <- brm(
  w1 ~ session * environment + (session * environment | user_id),
  data = params_happiness_clean_df,
  family = asym_laplace(),
  algorithm = "meanfield"
)
pp_check(m1) + xlim(-3, 3)
summary(m1)
conditional_effects(m1, "session")

m2 <- brm(
  w2 ~ session * environment + (session * environment | user_id),
  data = params_happiness_clean_df,
  family = asym_laplace(),
  algorithm = "meanfield"
)
pp_check(m2) + xlim(-3, 3)
summary(m2)
conditional_effects(m2, "session:environment")


m3 <- brm(
  w3 ~ session * environment + (session * environment | user_id),
  data = params_happiness_clean_df,
  family = student(),
  algorithm = "meanfield"
)
pp_check(m3) + xlim(-3, 3)
summary(m3)
conditional_effects(m3, "session")


m4 <- brm(
  w4 ~ session * environment + (session * environment | user_id),
  data = params_happiness_clean_df,
  family = student(),
  algorithm = "meanfield"
)
pp_check(m4) + xlim(-3, 3)
summary(m4)
conditional_effects(m4, "environment")


m5 <- brm(
  w5 ~ session * environment + (session * environment | user_id),
  data = params_happiness_clean_df,
  family = student(),
  algorithm = "meanfield"
)
pp_check(m5) + xlim(-3, 3)
summary(m5)
conditional_effects(m5, "environment")






