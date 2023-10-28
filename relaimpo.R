# Simulate some data
set.seed(123)
n <- 100
data <- data.frame(
  outcome = runif(n, 0, 1),
  stimulus = runif(n, -1, 1),
  RPE = runif(n, -1, 1),
  zmoodpre = runif(n, -1, 1),
  zcontrol = runif(n, -1, 1)
)

# Your initial parameters
params <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
names(params) <- c("w0", "w1", "w2", "w3", "w4", "w5", "gamma")  # Name the parameters

# Compute baseline happiness
baseline_happiness <- compute_predicted_happiness(params, data)

# Initialize a data frame to store the results
sensitivity_df <- data.frame()

# Loop through each parameter
for (i in 1:6) {
  # Small change in the parameter
  delta <- 0.01
  
  # Update one parameter at a time
  new_params <- params
  new_params[i] <- new_params[i] + delta
  
  # Compute new happiness
  new_happiness <- compute_predicted_happiness(new_params, data)
  
  # Compute the change in happiness
  delta_happiness <- mean(abs(new_happiness - baseline_happiness))
  
  # Store the results
  sensitivity_df <- rbind(
    sensitivity_df,
    data.frame(
      parameter = names(params)[i],
      delta_happiness = delta_happiness
    )
  )
}

# Sort by the change in happiness to see which parameter has the most impact
sensitivity_df <- sensitivity_df %>% arrange(desc(delta_happiness))






tar_load(params_happiness_clean_df)


# Add these needed variables
final_df$zmoodpre <- prl_df$zmoodpre
final_df$zcontrol <- prl_df$zcontrol



# Initialize a data frame to store the results for all user_ids
all_sensitivity_df <- data.frame()

# Get unique user_ids
unique_users <- unique(final_df$user_id)

# Loop through each unique user_id
for (user in unique_users) {
  # Extract data for this user
  real_data <- final_df %>%
    filter(user_id == user) %>%
    select(outcome, stimulus, rpe, happiness, zmoodpre, zcontrol)
  
  # Extract actual fitted parameters for this user
  # Assuming you have a data frame 'params_df' containing fitted params for each user
  actual_params <- params_happiness_clean_df %>%
    filter(user_id == user) %>%
    select("w0", "w1", "w2", "w3", "w4", "w5", "gamma", "user_id") |> 
    select(-user_id) %>% 
    unlist() %>% 
    as.numeric()
  
  names(actual_params) <- c("w0", "w1", "w2", "w3", "w4", "w5", "gamma")
  
  # Compute baseline happiness with real data and actual params
  baseline_happiness <- compute_predicted_happiness(actual_params, real_data)
  
  # Initialize a data frame to store the results for this user
  sensitivity_df <- data.frame()
  
  # Loop through each parameter
  for (i in 1:6) {
    # Small change in the parameter
    delta <- 0.01
    
    # Update one parameter at a time
    new_params <- actual_params
    new_params[i] <- new_params[i] + delta
    
    # Compute new happiness with updated params
    new_happiness <- compute_predicted_happiness(new_params, real_data)
    
    # Compute the change in happiness
    delta_happiness <- mean(abs(new_happiness - baseline_happiness))
    
    # Store the results
    sensitivity_df <- rbind(
      sensitivity_df,
      data.frame(
        user_id = user,  # Add user_id for tracking
        parameter = names(actual_params)[i],
        delta_happiness = delta_happiness
      )
    )
  }
  
  # Combine with the results for other users
  all_sensitivity_df <- rbind(all_sensitivity_df, sensitivity_df)
}

# Sort by the change in happiness to see which parameter has the most impact
sorted_sensitivity_df <- all_sensitivity_df %>% arrange(desc(delta_happiness))

filtered_df <- sorted_sensitivity_df %>%
  filter(if_any(everything(), ~ . <= 3))

filtered_df |> 
  ggplot(aes(x = parameter, y = delta_happiness)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    theme_minimal() +
    ggtitle('Sensitivity Analysis of Happiness Model Parameters') +
    xlab('Parameters') +
    ylab('Change in Predicted Happiness')

df$parameter <- as.factor(df$parameter)

rf_model <- randomForest(
  delta_happiness ~ parameter, data = df, ntree = 100)

importance(rf_model)
varImpPlot(rf_model)




# Calculate the IQR for delta_happiness
Q1 <- quantile(sorted_sensitivity_df$delta_happiness, 0.25)
Q3 <- quantile(sorted_sensitivity_df$delta_happiness, 0.75)
IQR <- Q3 - Q1

# Define bounds for the outliers
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Remove outliers
filtered_df <- subset(sorted_sensitivity_df, delta_happiness >= lower_bound & delta_happiness <= upper_bound)


wide_df <- filtered_df %>%
  pivot_wider(names_from = parameter, values_from = delta_happiness)


full_model <- lm(delta_happiness ~ parameter, data = filtered_df)
reduced_model1 <- lm(delta_happiness ~ parameter - w0, data = filtered_df)

data_without_w1 <- subset(filtered_df, parameter != "w1")
reduced_model1 <- lm(delta_happiness ~ parameter, data = data_without_w1)

data_without_w1_w2 <- subset(filtered_df, !(parameter %in% c("w1", "w2")))
reduced_model2 <- lm(delta_happiness ~ parameter, data = data_without_w1_w2)

data_without_w1_w2_w3 <- subset(filtered_df, !(parameter %in% c("w1", "w2")))
reduced_model2 <- lm(delta_happiness ~ parameter, data = data_without_w1_w2)

aic_full_model <- AIC(full_model)
bic_full_model <- BIC(full_model)

aic_full_model <- AIC(reduced_model1)
bic_full_model <- BIC(reduced_model1)

# Fit the model
fit <- brm(
  formula = delta_happiness ~ 1 + (1|user_id) + (1|parameter), 
  data = filtered_df, 
  family = asym_laplace(),
  prior = c(
    set_prior("normal(0, 5)", class = "Intercept"),
    set_prior("normal(0, 5)", class = "sd"),
    set_prior("normal(0, 5)", class = "sigma")
  ),
  iter = 2000, 
  algorithm = "meanfield",
  chains = 4
)

filtered_df$dh <- as.vector(scale(filtered_df$delta_happiness))

fit <- brm(
  dh ~ parameter + (parameter | user_id), 
  family = student(),
  backend = "cmdstanr",
  data = filtered_df
)
pp_check(fit) + xlim(-5, 5)

random_effects <- ranef(fit)$parameter[, , "Estimate"]

summary_stats <- data.frame(
  Mean = rowMeans(random_effects),
  SD = apply(random_effects, 1, sd)
)

ggplot(data = summary_stats, aes(x = rownames(summary_stats), y = Mean)) +
  geom_point(aes(size = SD)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  ggtitle("Random Effects for Parameters") +
  xlab("Parameter") +
  ylab("Posterior Mean of Random Effects")




