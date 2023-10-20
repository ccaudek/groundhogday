
# get parameters

# all subjects and sessions parameters
params <- results_df

# Select one subject
id <- user_id_codes[50]


# da qui!!
onesubj_data <- dz_clean |> 
  dplyr::filter(user_id == id)

n_ema_episodes <- length(unique(onesubj_data$ema_number))

i <- n_ema_episodes

par_list <- list()

for (i in seq_len(n_ema_episodes)) {
  
  ema_session <- onesubj_data |> 
    dplyr::filter(ema_number == i)
  
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
  
  subj_code <- unique(ema_session$user_id)
  
  subj_session_params <- params |> 
    dplyr::filter(user_id == subj_code & ema_number == i) |> 
    dplyr::select("w0", "w1", "w2", "w3", "w4", "w5", "gamma")
  
  happiness_hat <- compute_predicted_happiness(
    subj_session_params, df
  )
  
  par_list[[i]] <- happiness_hat
}

# TODO
# I have a list of ema_number elements. Each element is a vector of 30 values,
# which correspond to the predicted momentary happiness in each trial of that
# session.

# Convert list to dataframe
par_df <- do.call(rbind, par_list)
par_df <- as.data.frame(par_df)
colnames(par_df) <- c(
  "w0", "w1", "w2", "w3", "w4", "gamma", 
  "is_reversal", "ema_number", "user_id"
)

plot(df$happiness, happiness_hat)


plot(1:30, df$happiness, type = "l", col = "black")
lines(1:30, happiness_hat, type = "l", col = "blue")

cor(df$happiness, happiness_hat, method = "spearman")




################################################################################

# BELOW
# ChatGPT scripts for the RPE

# REQUIRES the df DataFrame for a subject and session:
#
# $ trial               <int> 1, 2, 3, 4, 5, 
# $ stimulus            <dbl> -1, 1, 1, 1, 1, 
# $ reversal            <dbl> 0, 0, 0, 0, 0, 
# $ outcome             <dbl> -1, 1, 1, 1, 
# $ delta_p             <dbl> 0.0, 0.0, 0.0, 
# $ happiness           <dbl> -0.8576963, -2.4460228, 






# Compute alpha ----------------------------------------------------------------

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


add_pre <- function(df, best_alpha) {
  
  # Assuming best_alpha contains the estimated alpha value
  # best_alpha <- result$par  # Uncomment this line if you're running this 
  # after the MLE code
  
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

# Compute RPE ------------------------------------------------------------------

# Assuming best_alpha contains the estimated alpha value
# best_alpha <- result$par  # Uncomment this line if you're running this 
# after the MLE code

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

# Print out the updated data frame
print(df)











