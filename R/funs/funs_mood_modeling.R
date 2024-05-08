#' Groundhog Day project.
#'
#' Functions to estimate the parameters of the momentary happiness model 
#' for the data of the PRL task of the Groundhog Day project.


#' get_data() ------------------------------------------------------------------
#' @description
#' Read pre-processed complete raw PRL data and remove bad subjects. Do some
#' minimial housekeeping (standardize by subject the critical variables).
#' The data have already been cleaned so that each participant completed
#' at least 4 EMA sessions.
#' @param file_path The path to the pre-processed RDS file with all raw data.
#' @returns A data frame.
#' 
get_data <- function(file_path) {
  # List of bad user_ids
  bad_ids <- c(
    3338029881, 3665345709, 3248648540, 3668615349, 3456996255,
    3475648095, 3497824506, 3888185906, 3667000898
  )
  
  # Read, filter, and transform the data in a single chain
  data <- readRDS(file_path) %>%
    dplyr::mutate(
      date = if ("date" %in% names(.)) as.Date(date) else date,
      user_id = as.numeric(user_id)
    ) %>%
    dplyr::filter(
      !is.na(user_id),
      !(user_id %in% bad_ids),
      ema_number < 13
    ) %>%
    dplyr::group_by(user_id) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      ema_number = dense_rank(date),
      zmoodpre = as.vector(scale(mood_pre, center = TRUE, scale = TRUE)),
      zcontrol = as.vector(scale(control, center = TRUE, scale = TRUE)),
      zim = as.vector(scale(instant_mood, center = TRUE, scale = TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      zcontrol = ifelse(is.na(zcontrol), 0, zcontrol)
    ) %>%
    dplyr::select(
      instant_mood, feedback, trial, user_id, mood_pre, control, mood_post,
      is_reversal, is_target_chosen, ema_number, zim, zmoodpre, zcontrol,
      date
    )
  
  return(data)
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
    group_by(user_id, ema_number, trial) %>%  # Group by trial as well
    summarize(
      mood_pre = mean(mood_pre),
      mood_post = mean(mood_post),
      control = mean(control),
      .groups = 'drop'  # Drop the grouping
    )
  
  # Final data frame.
  results_df <- left_join(
    all_results_df, bysubj_mood_df,
    by = c("user_id", "ema_number", "trial")  # Join by trial as well
  )
  
  return(results_df)
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
  result_dfs <- list()
  
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
    init_params <- c(0, 0, 0, 0, 0, 0, 0, 0.5)
    
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
    
    # Compute the predicted happiness values
    predicted_happiness <- compute_predicted_happiness(mle_params, df)
    
    # Add the estimated parameters and predicted happiness to each row
    df$predicted_happiness <- predicted_happiness
    df$mle_params <- list(mle_params)
    
    # Add session-level variables to each row
    df$is_reversal <- ifelse(unique(ema_session$is_reversal) == "yes", 1, 0)
    df$ema_number <- i
    df$user_id <- id
    df$best_alpha <- best_alpha
    
    return(df)
  }
  
  # Loop through each EMA session
  for (i in seq_len(n_ema_episodes)) {
    result_dfs[[i]] <- process_ema_session(i)
  }
  
  # Bind all data frames together into a single data frame
  result_df <- do.call(rbind, result_dfs)
  
  return(result_df)
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
  w6 <- params[7] # weight for trial
  gamma <- params[8] # weight for gamma

  zmoodpre <- data$zmoodpre
  zcontrol <- data$zcontrol
  trial <- data$trial

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
      w5 * zcontrol[t] +
      w6 * trial[t] 
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

#' detect_outliers() -----------------------------------------------------------
#' @description
#' Helper function to detect outliers based on Mahalanobis distance
#' 
detect_outliers <- function(df, cols, prob = 0.99) {
  if (ncol(df) < 1 || nrow(df) < 1) {
    stop("Data frame must have at least one row and one column.")
  }
  
  set.seed(123)
  values <- df %>% dplyr::select(all_of(cols))
  
  center <- colMeans(values)
  cov_matrix <- cov(values)
  
  values$mdist <- mahalanobis(x = values, center = center, cov = cov_matrix)
  
  cutoff <- qchisq(p = prob, df = ncol(values))
  
  values %>%
    mutate(is_outlier = ifelse(mdist > cutoff, 1, 0))
}

#' perform_imputation() --------------------------------------------------------
#' @description
#' Helper function to perform multiple imputation
#' 
perform_imputation <- function(df, cols, m = 1, seed = 500) {
  imputed_data <- mice(
    df %>% dplyr::select(all_of(cols)), 
    m = m, 
    maxit = 50, 
    method = "pmm", 
    seed = seed
  )
  complete(imputed_data, 1)
}

#' clean_params_happiness_model() ----------------------------------------------
#' @description
#' Main function to clean parameters related to the happiness model
#'
clean_params_happiness_model <- function(params_happiness_df) {
  
  # Helper function to rename MLE parameters
  rename_mle_params <- function(df) {
    df %>% 
      rename(
        w0 = mle_params_1,
        w_outcome = mle_params_2,
        w_stimulus = mle_params_3,
        w_rpe = mle_params_4,
        w_moodpre = mle_params_5,
        w_control = mle_params_6,
        w_trial = mle_params_7,
        w_gamma = mle_params_8
      )
  }
  
  # Step 1: Unnest and rename MLE parameters
  expanded_df <- params_happiness_df %>%
    unnest_wider(mle_params, names_sep = "_") %>%
    rename_mle_params()
  
  # Step 1.1: Impute outliers for MLE parameters
  outlier_cols1 <- c("w0", "w_outcome", "w_stimulus", "w_rpe", "w_moodpre", "w_control", "w_trial", "w_gamma")
  imputed_df1 <- impute_outliers_params_happiness(expanded_df, outlier_cols1)
  
  # Step 2: Impute outliers for mood and control variables
  outlier_cols2 <- c("mood_pre", "control", "mood_post")
  imputed_df2 <- impute_outliers_mood(imputed_df1, outlier_cols2)
  
  # Step 3: Additional transformations
  final_df <- imputed_df2 %>%
    group_by(user_id) %>%
    mutate(
      mood_dif = mood_post - mood_pre,
      mood_pre_cw = mood_pre - mean(mood_pre, trim = 0.1),
      environment = ifelse(is_reversal == 1, "Volatile", "Stable")
    ) %>%
    ungroup() |> 
    dplyr::rename(alpha = best_alpha)
  
  return(final_df)
}


#' impute_outliers_params_happiness() ------------------------------------------
#' @description
#' Function to impute outliers for happiness parameters
#' 
impute_outliers_params_happiness <- function(df, cols) {
  outliers <- detect_outliers(df, cols)
  outlier_indices <- which(outliers$is_outlier == 1)
  df[outlier_indices, cols] <- NA
  imputed_data <- perform_imputation(df, cols)
  df[cols] <- imputed_data
  return(df)
}

#' impute_outliers_mood() ------------------------------------------------------
#' @description
#' Function to impute outliers for mood and control variables
#' 
impute_outliers_mood <- function(df, cols) {
  outliers <- detect_outliers(df, cols)
  outlier_indices <- which(outliers$is_outlier == 1)
  df[outlier_indices, cols] <- NA
  imputed_data <- perform_imputation(df, cols)
  df[cols] <- imputed_data
  return(df)
}


# eof ----

