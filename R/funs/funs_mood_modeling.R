# Functions for estimating the parameters of the momentary happiness model 
# for the data of the PRL task of the Groundhog Day project.


#' remove_outliers() -----------------------------------------------------------
#' 
#' @description
#' Remove with NAs the values above a threshold.
#' @param my_vector A vector.
#' @return A vector.
#' 
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
#' Read pre-processed complete raw data and remove bad subjects. Do some
#' minimial housekeeping (standardize by subject the critical variables).
#' The data have already been cleaned so that each participant completed
#' at least 4 EMA sessions.
#' @param file_path The path to the pre-processed RDS file with all raw data.
#' @returns A data frame.
#' 
get_data <- function(file_path) {

  # file_path <- here::here("data", "prep", "groundhog_all_clean.RDS")
  d1 <- readRDS(file_path)

  if ("date" %in% names(d1)) {
    d1$date <- as.Date(d1$date)
  }

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

  # table(d3$ema_number)
  #    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
  # 6450 6450 6450 6450 6090 5820 5550 5250 4890 4170 3600 2760 1890  390   60   30   30

  # Remove the last sessions having only few participants.
  d4 <- d3 |>
    dplyr::filter(ema_number < 13)

  # Compute ema_number (rank-order of session for each participant) and
  # standardize instant_mood, mood_pre, mood_post, and control by user_id.
  d5 <- d4 |>
    group_by(user_id) |>
    arrange(date) |>
    mutate(
      ema_number = dense_rank(date),
      zmoodpre = as.vector(scale(mood_pre, center = TRUE, scale = TRUE)),
      zcontrol = as.vector(scale(control, center = TRUE, scale = TRUE)),
      zim = as.vector(scale(instant_mood, center = TRUE, scale = TRUE))
    ) |>
    ungroup()

  # For some user_id, control has the same value in all EMA sessions and,
  # therefore, zcontrol is NA. Replace these NAs with 0.
  d5$zcontrol <- ifelse(is.na(d5$zcontrol), 0, d5$zcontrol)

  # # Discard the first session (i.e., ema_number == 1) for each participant.
  # d5_filtered <- d5 %>%
  #   dplyr::filter(ema_number != 1)
  # 
  # # Recompute ema_number starting from 1 for each participant.
  # d5_recomputed <- d5_filtered %>%
  #   dplyr::group_by(user_id) %>%
  #   dplyr::arrange(date) %>%
  #   dplyr::mutate(
  #     ema_number = dense_rank(date)
  #   ) %>%
  #   dplyr::ungroup()

  # Select the variables used for the momentary mood modeling.
  d6 <- d5 |>
    dplyr::select(
      instant_mood, feedback, trial, user_id, mood_pre, control, mood_post,
      is_reversal, is_target_chosen, ema_number, zim, zmoodpre, zcontrol,
      date
    )

  return(d6)
}


# get_params_happiness_model() -------------------------------------------------

#' @description
#' Estimate the parameters of the momentary happiness model for each 
#' participant and each EMA session. 
#' 
#' @param dz_clean A data frame.
#' @param user_id_codes A vector.
#' @returns A data frame with the estimated parameters, together with 
#' control, mood_pre, and mood_post.
#' 
get_params_happiness_model <- function(df, user_id_codes) {
  
  # Apply the function process_user() to each user_id.
  results_list <- lapply(user_id_codes, process_user, df)
  
  # Bind all data frames together into a single data frame.
  all_results_df <- bind_rows(results_list)
  
  # Add further needed variables.
  bysubj_mood_df <- df %>%
    group_by(user_id, ema_number) %>%
    summarize(
      mood_pre = mean(mood_pre),
      mood_post = mean(mood_post),
      control = mean(control),
      zmoodpre = mean(zmoodpre),
      zcontrol = mean(zcontrol)
    ) %>%
    ungroup()
  
  # Final data frame.
  results_df <- left_join(
    all_results_df, bysubj_mood_df,
    by = c("user_id", "ema_number")
  )
  
  results_df
}


#' process_user() --------------------------------------------------------------
#' 
#' @param id The string `user-id` for a single participant. 
#' @param dat A data frame with the complete data set.
#' 
process_user <- function(id, dat) {
  # Set random seed for reproducibility
  set.seed(12345)
  
  # Filter data for the current user_id
  onesubj_data <- dat %>% dplyr::filter(user_id == id)
  
  # Compute best_alpha once for each user_id
  best_alpha <- get_alpha_all_sessions(onesubj_data)
  
  # Count the number of unique EMA episodes for this user
  n_ema_episodes <- length(unique(onesubj_data$ema_number))
  
  # Initialize the list to hold results
  par_list <- list()
  
  # Function to process a single EMA session
  process_ema_session <- function(i) {
    ema_session <- onesubj_data %>%
      dplyr::filter(ema_number == i) %>%
      dplyr::select(
        user_id, ema_number, trial, is_target_chosen, is_reversal,
        feedback, zim, zmoodpre, zcontrol
      )
    
    # Create a new data frame with all the required information for one single session
    df <- data.frame(
      trial = ema_session$trial,
      stimulus = ifelse(ema_session$is_target_chosen == 0, -1, ema_session$is_target_chosen),
      reversal = ifelse(ema_session$is_reversal == "yes", c(rep(0, 15), 1, rep(0, 14)), rep(0, 30)),
      outcome = ifelse(ema_session$feedback == 0, -1, ema_session$feedback),
      delta_p = ifelse(ema_session$is_reversal == "yes", c(rep(0, 15), 0.6, rep(0, 14)), rep(0, 30)),
      happiness = ema_session$zim,
      zmoodpre = ema_session$zmoodpre,
      zcontrol = ema_session$zcontrol
    )
    
    # Add the RPE column using the pre-computed best_alpha
    df <- add_rpe(df, best_alpha)
    
    # Initial parameter estimates
    init_params <- c(0, 0, 0, 0, 0, 0, 0.5)
    
    # Optimization with try-catch
    opt_result <- tryCatch(
      {
        optim(init_params, nll, data = df)
      },
      error = function(e) {
        warning("Optimization failed for this set of parameters, returning NA")
        return(NULL)
      }
    )
    
    # Extract the parameters or set them to NA if optimization failed
    if (is.null(opt_result)) {
      mle_params <- rep(NA, length(init_params))
    } else {
      mle_params <- opt_result$par
    }
    
    # Create the output vector
    out <- c(
      mle_params,
      ifelse(unique(ema_session$is_reversal) == "yes", 1, 0),
      i,
      id,
      best_alpha
    )
    
    return(out)
  }
  
  # Parallelize the loop using mclapply
  par_list <- mclapply(seq_len(n_ema_episodes), process_ema_session, mc.cores = 4)
  
  # Convert the list to a data frame
  par_df <- do.call(rbind, par_list)
  par_df <- as.data.frame(par_df)
  colnames(par_df) <- c(
    "w0", "w1", "w2", "w3", "w4", "w5", "gamma", 
    "is_reversal", "ema_number", "user_id", "alpha"
  )
  
  return(par_df)
}


#' get_alpha() -----------------------------------------------------------------
#' 
#' @description
#' Compute alpha for a single session of the current user_id.
#' @param df the data frame for a single EMA session of the current user_id. 
#' @return alpha
#' 
get_alpha <- function(df) {
  neg_log_likelihood <- function(alpha, df) {
    predictions <- list("1" = 0, "-1" = 0)
    nll <- 0

    for (i in 1:nrow(df)) {
      row <- df[i, ]
      stimulus <- as.character(row$stimulus)
      outcome <- row$outcome

      prediction <- predictions[[stimulus]]

      # Calculate RPE
      rpe <- outcome - prediction

      # Update prediction for the stimulus
      new_prediction <- prediction + alpha * rpe

      # Store the new prediction back into the list
      predictions[[stimulus]] <- new_prediction

      # Accumulate the negative log-likelihood (using Gaussian distribution for simplicity)
      nll <- nll - log(dnorm(outcome, mean = prediction, sd = 1))
    }

    return(nll)
  }

  result <- optim(
    par = 0.01, # initial value for alpha
    fn = neg_log_likelihood,
    df = df,
    method = "Brent",
    lower = 0,
    upper = 1
  )

  best_alpha <- result$par
  best_alpha
}


#' get_alpha_all_sessions() ----------------------------------------------------
#' 
#' @description
#' Compute alpha for all session of the current user_id.
#' @param df the data frame for a single EMA session of the current user_id. 
#' @return alpha
#' 
get_alpha_all_sessions <- function(onesubj_data) {
  
  # Create a new data frame with all the required information for the current 
  # participant.
  df <- data.frame(
    trial = onesubj_data$trial,
    # 1 for stimulus "same action", -1 for stimulus "different action"
    stimulus = ifelse(
      onesubj_data$is_target_chosen == 0, -1, onesubj_data$is_target_chosen
    ), 
    outcome = ifelse(onesubj_data$feedback == 0, -1, onesubj_data$feedback)
  )
  
  neg_log_likelihood <- function(alpha, df) {
    predictions <- list("1" = 0, "-1" = 0)
    nll <- 0
    
    for (i in 1:nrow(df)) {
      row <- df[i, ]
      stimulus <- as.character(row$stimulus)
      outcome <- row$outcome
      
      prediction <- predictions[[stimulus]]
      
      # Calculate RPE
      rpe <- outcome - prediction
      
      # Update prediction for the stimulus
      new_prediction <- prediction + alpha * rpe
      
      # Store the new prediction back into the list
      predictions[[stimulus]] <- new_prediction
      
      # Accumulate the negative log-likelihood (using Gaussian distribution for simplicity)
      nll <- nll - log(dnorm(outcome, mean = prediction, sd = 1))
    }
    
    return(nll)
  }
  
  result <- optim(
    par = 0.01, # initial value for alpha
    fn = neg_log_likelihood,
    df = df,
    method = "Brent",
    lower = 0,
    upper = 1
  )
  
  best_alpha <- result$par
  best_alpha
}


#' add_rpe() -------------------------------------------------------------------
#' 
#' @description
#' Adds the Reward Prediction Error to the df data frame for the single session 
#' of the current user_id. 
#' @param df A data frame.
#' @param best_alpha The estimated alpha for the current EMA session.
#' @return the df data frame with added the RPE column.
#' 
add_rpe <- function(df, best_alpha) {
  # Initialize predictions for each stimulus type
  predictions <- list("1" = 0, "-1" = 0)

  # Create empty vectors to store the RPE and updated predictions
  rpe_vector <- numeric()
  updated_predictions_vector <- numeric()

  # Loop through each row of the data frame
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    stimulus <- as.character(row$stimulus)
    outcome <- row$outcome

    # Retrieve the current prediction for the stimulus
    prediction <- predictions[[stimulus]]

    # Calculate the RPE
    rpe <- outcome - prediction
    rpe_vector <- c(rpe_vector, rpe)

    # Update the prediction for the stimulus
    new_prediction <- prediction + best_alpha * rpe
    updated_predictions_vector <- c(updated_predictions_vector, new_prediction)

    # Update the predictions list
    predictions[[stimulus]] <- new_prediction
  }

  # Add the RPE and updated predictions to the data frame
  df$RPE <- rpe_vector
  df$Updated_Predictions <- updated_predictions_vector

  df
}


# compute_predicted_happiness() ------------------------------------------------
#' @description
#' Compute the predicted happiness given parameters and data.
#'
#' @param params A vector [w0, w1, w2, w3, w4, w5, gamma].
#' @param data A dataframe with columns: outcome, reversal, stimulus, delta_p,
#' RPE.
#' @return A vector of 30 elements.
#' @examples
#' compute_predicted_happiness(params, df)
#'
compute_predicted_happiness <- function(params, data, ...) {
  w0 <- params[1] # intercept
  w1 <- params[2] # weight for outcome
  w2 <- params[3] # weight for the chosen stimulus
  w3 <- params[4] # weight for RPE
  w4 <- params[5] # weight for zmoodpre
  w5 <- params[6] # weight for zcontrol
  gamma <- params[7] # weight for gamma

  zmoodpre <- data$zmoodpre
  zcontrol <- data$zcontrol

  predicted_happiness <- numeric(nrow(data))

  for (t in 1:nrow(data)) {
    weighted_outcome <- sum(gamma^(0:(t - 1)) * rev(data$outcome[1:t]))
    weighted_stimulus <- sum(gamma^(0:(t - 1)) * rev(data$stimulus[1:t]))
    weighted_RPE <- sum(gamma^(0:(t - 1)) * rev(data$RPE[1:t]))

    predicted_happiness[t] <- w0 +
      w1 * weighted_outcome +
      w2 * weighted_stimulus +
      w3 * weighted_RPE +
      w4 * zmoodpre[t] +
      w5 * zcontrol[t]
  }

  return(predicted_happiness)
}


# nll() ------------------------------------------------------------------------

#' @description
#' Compute the negative log-likelihood assuming normally distributed residuals.
#' Function to be used with optim().
#' 
#' @param params A vector [w0, w1, w2, w3, w4, w5, gamma].
#' @param data A dataframe.
#' 
nll <- function(params, data) {
  predicted_happiness <- compute_predicted_happiness(params, data)
  ssr <- sum((data$happiness - predicted_happiness)^2)
  n <- nrow(data)
  nll_value <- n/2 * log(2 * pi) + n/2 * log(ssr/n) + n/2
  return(nll_value)
}


#' clean_params_happiness_model() ----------------------------------------------
#' 
#' @description
#' Find outiers in the parameters of the momentary happiness model, replace
#' them with NAs and perform multiple imputation.
#' @param params_happiness_df A data frame.
#' @return A data frame.
#' 
clean_params_happiness_model <- function(params_happiness_df) {
  
  # Outlier detection and multiple imputation of the parameters of the
  # momentary happiness model's parameters.
  df <- impute_outliers_params_happiness(params_happiness_df)
  # Outlier detection and multiple imputation on the mood and control
  # variables.
  df1 <- impute_outliers_mood(df)

  # On the imputed data, compute mood_dif and mood_pre_cw.
  df2 <- df1 |>
    group_by(user_id) |>
    mutate(
      mood_dif = mood_post - mood_pre,
      mood_pre_cw = mood_pre - mean(mood_pre, trim = 0.1)
    ) |> 
    ungroup()

  df2$environment <- ifelse(df2$is_reversal == 1, "Volatile", "Stable")

  return(df2)
}


#' impute_outliers_params_happiness() ------------------------------------------
#' 
#' @description
#' Use Mahalanobis distance to identify outliers in the estimates of the 
#' momentary happiness model's parameters. Replace outliers with NAs. Perform
#' multiple imputation on the missing values. 
#' @param params_happiness_df A data frame.
#' @return The imputed params_happiness_df data frame.
#' 
impute_outliers_params_happiness <- function(params_happiness_df) {
  PROB <- 0.99
  IMPUTATION_M <- 5
  
  set.seed(123)
  
  # Step 1: Detect Outliers
  params_values <- params_happiness_df %>%
    dplyr::select(w0, w1, w2, w3, w4, w5, gamma)

  params_center <- colMeans(params_values)
  params_cov <- cov(params_values)

  params_values$mdist <- mahalanobis(
    x = params_values,
    center = params_center,
    cov = params_cov
  )

  cutoff <- qchisq(p = PROB, df = ncol(params_values[, 1:7]))

  params_values <- params_values %>%
    mutate(is_outlier = ifelse(mdist > cutoff, 1, 0))

  outlier_indices <- which(params_values$is_outlier == 1)

  # Replace outliers with NA
  params_happiness_df[
    outlier_indices, c("w0", "w1", "w2", "w3", "w4", "w5", "gamma")
  ] <- NA

  # Step 2: Multiple Imputation
  # Select only the numeric variables of interest
  temp <- params_happiness_df |>
    dplyr::select(
      w0, w1, w2, w3, w4, w5, gamma, is_reversal, ema_number, alpha,
      zmoodpre, zcontrol
    )
  imputed_data <- mice(
    temp, m = IMPUTATION_M, maxit = 50, method = "pmm", seed = 500
  )
  completed_data <- complete(imputed_data, 1) # using the first imputed dataset

  params_happiness_df[, names(completed_data)] <- completed_data

  return(params_happiness_df)
}


#' impute_outliers_mood() ------------------------------------------------------
#' 
#' @description
#' Outlier detection and multiple imputation for the mood and control variables.
#' @param params_happiness_df A data frame.
#' @return A data frame.
#' 
impute_outliers_mood <- function(params_happiness_df) {
  set.seed(123)
  # Step 1: Detect Outliers
  mood_values <- params_happiness_df %>%
    dplyr::select(mood_pre, control, mood_post)
  
  mood_center <- colMeans(mood_values)
  mood_cov <- cov(mood_values)
  
  mood_values$mdist <- mahalanobis(
    x = mood_values,
    center = mood_center,
    cov = mood_cov
  )
  
  cutoff <- qchisq(p = 0.99, df = ncol(mood_values[, 1:3]))
  
  mood_values <- mood_values %>%
    mutate(is_outlier = ifelse(mdist > cutoff, 1, 0))
  
  outlier_indices <- which(mood_values$is_outlier == 1)
  
  # Replace outliers with NA
  params_happiness_df[
    outlier_indices, c("mood_pre", "control", "mood_post")
  ] <- NA
  
  # Step 2: Multiple Imputation
  # Select only the numeric variables of interest
  temp <- params_happiness_df |>
    dplyr::select(mood_pre, control, mood_post)
  imputed_data <- mice(temp, m = 5, maxit = 50, method = "pmm", seed = 500)
  completed_data <- complete(imputed_data, 1) # using the first imputed dataset
  
  params_happiness_df[, names(completed_data)] <- completed_data
  
  return(params_happiness_df)
}


# eof ----

