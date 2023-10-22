# Functions for computing the parameters of the model accounting for
# momentary happiness in a PRL task.


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
  w0 <- params[1]
  w1 <- params[2]
  w2 <- params[3]
  w3 <- params[4]
  w4 <- params[5]
  w5 <- params[6]
  gamma <- params[7]
  
  predicted_happiness <- numeric(nrow(data))
  
  for (t in 1:nrow(data)) {
    predicted_happiness[t] <- w0 +
      w1 * sum(gamma^(t-1:(t-1)) * data$outcome[1:(t-1)]) +
      w2 * sum(gamma^(t-1:(t-1)) * data$reversal[1:(t-1)]) +
      w3 * sum(gamma^(t-1:(t-1)) * data$stimulus[1:(t-1)]) +
      w4 * sum(gamma^(t-1:(t-1)) * data$delta_p[1:(t-1)]) +
      w5 * sum(gamma^(t-1:(t-1)) * data$RPE[1:(t-1)])
  }
  
  res <- unlist(predicted_happiness)
  
  return(res)
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


# process_user() ---------------------------------------------------------------

# Define function to process a single user
process_user <- function(id, dz_clean) {
  
  set.seed(123)
  
  onesubj_data <- dz_clean |> 
    dplyr::filter(user_id == id)
  
  n_ema_episodes <- length(unique(onesubj_data$ema_number))
  
  par_list <- list()
  
  for (i in seq_len(n_ema_episodes)) {
    
    ema_session <- onesubj_data |> 
      dplyr::filter(ema_number == i) |> 
      dplyr::select(
        user_id, ema_number, trial, is_target_chosen, is_reversal, 
        feedback, zim
      )
    
    # Required information for a single session of a subject.
    df <- data.frame(
      trial = ema_session$trial,
      stimulus = ifelse(
        ema_session$is_target_chosen == 0, -1, ema_session$is_target_chosen
      ), # 1 for stimulus A, -1 for stimulus B
      # reversal = c(rep(0, 14), 1, rep(0, 15)),  # Reversal occurs at trial 15
      reversal = ifelse(
        ema_session$is_reversal == "yes",
        c(rep(0, 15), 1, rep(0, 14)),
        rep(0, 30)
      ),
      outcome = ifelse(ema_session$feedback == 0, -1, ema_session$feedback),
      # delta_p = c(rep(0, 15), 0.6, rep(0, 14)),  # Change in probability at reversal
      delta_p = ifelse(
        ema_session$is_reversal == "yes",
        c(rep(0, 15), 0.6, rep(0, 14)),
        rep(0, 30)
      ),
      happiness = ema_session$zim # standardized by user_id
    )
    
    # Get alpha for a single ema_number session and the current user_id
    best_alpha = get_alpha(df)
    # Add the RPE column to the df DataFrame
    df = add_pre(df, best_alpha)
    
    # Optimize
    # Initial guesses for w0, w1, w2, w3, w4, w5, and gamma
    init_params <- c(0, 0, 0, 0, 0, 0, 0.5)  
    opt_result <- optim(init_params, nll, data=df)
    mle_params <- opt_result$par
    # add further information
    out <- c(
      mle_params, 
      ifelse(unique(ema_session$is_reversal) == "yes", 1, 0), # is_reversal
      i, # ema_number
      id, # user_id
      best_alpha # alpha
    )
    
    par_list[[i]] <- out
  }
  
  # Convert list to DataFrame
  par_df <- do.call(rbind, par_list)
  par_df <- as.data.frame(par_df)
  colnames(par_df) <- c(
    "w0", "w1", "w2", "w3", "w4", "w5", "gamma", 
    "is_reversal", "ema_number", "user_id", "alpha"
  )
  
  # cat('user_id:', unique(onesubj_data$user_id), '\n')
  
  return(par_df)
}


# get_alpha() ------------------------------------------------------------------

#' @description
#' Compute alpha for a single session of a subject.
#' @param df the DataFrame for a subject and session structured like this:
#'
#' $ trial               <int> 1, 2, 3, 4, 5, 
#' $ stimulus            <dbl> -1, 1, 1, 1, 1, 
#' $ reversal            <dbl> 0, 0, 0, 0, 0, 
#' $ outcome             <dbl> -1, 1, 1, 1, 
#' $ delta_p             <dbl> 0.0, 0.0, 0.0, 
#' $ happiness           <dbl> -0.8576963, -2.4460228, 
#' 
#' @return alpha

get_alpha <- function(df) {
  
  neg_log_likelihood <- function(alpha, df) {
    predictions <- list('1' = 0, '-1' = 0)
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
  
  result <- optim(par = 0.01,  # initial value for alpha
                  fn = neg_log_likelihood,
                  df = df,
                  method = "Brent",
                  lower = 0,
                  upper = 1)
  
  best_alpha <- result$par
  best_alpha
}


# add_pre() --------------------------------------------------------------------

#' @description
#' Adds the Reward Prediction Error to the df DataFrame for a single session 
#' of a subject.
#' @param df 
#' @param best_alpha
#' @return the df DataFrame with the RPE column

add_pre <- function(df, best_alpha) {
  
  # Initialize predictions for each stimulus type
  predictions <- list('1' = 0, '-1' = 0)
  
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


# eof ----



