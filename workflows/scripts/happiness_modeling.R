# Script for instant mood analysis
# Project: groundhog_day


# Initialize an empty data frame to store the results
final_df <- data.frame()

# Get all unique user_ids
unique_users <- unique(dz_clean$user_id)

# Loop through each unique user_id
for (id in unique_users) {
  onesubj_data <- dz_clean %>% dplyr::filter(user_id == id)
  n_ema_episodes <- length(unique(onesubj_data$ema_number))
  
  # Loop through each ema_number for the current user_id
  for (i in seq_len(n_ema_episodes)) {
    ema_session <- onesubj_data %>% dplyr::filter(ema_number == i)
    
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
    
    # After calculating happiness_hat, create a temporary data frame
    temp_df <- data.frame(
      user_id = id,
      ema_session = i,
      trial = df$trial,
      reversal = ema_session$is_reversal,
      stimulus = df$stimulus,
      outcome = df$outcome,
      rpe = df$RPE,
      happiness = df$happiness,
      happiness_hat = happiness_hat
    )
    
    # Append this data to the final data frame
    final_df <- rbind(final_df, temp_df)
  }
}


# Now final_df contains data for all users and all sessions

foo <- final_df |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1),
    h_hat = mean(happiness_hat, trim = 0.1)
  )
foo$reversal <- as.factor(foo$reversal)

foo |> 
  ggplot(aes(x=trial, y=h)) +
  geom_line() +
  facet_wrap(~ reversal)


final_df$h <- ifelse(
  final_df$happiness < -5 | final_df$happiness > 5, NA, 
  final_df$happiness
)

final_df$h_hat <- ifelse(
  final_df$happiness_hat < -15 | final_df$happiness_hat > 15, NA, 
  final_df$happiness_hat
)

cor(final_df$happiness, final_df$h_hat, use="pairwise.complete.obs")

# Correlation between momentary earning and momentary happiness
temp <- final_df
temp <- temp |> 
  group_by(user_id, ema_session) |> 
  mutate(
    earning = cumsum(outcome)
  )
cor(temp$happiness, temp$earning, use="pairwise.complete.obs")

fit <- brm(
  happiness ~ h_hat + earning + (1 + h_hat + earning | user_id/ema_session),
  data = temp,
  family = asym_laplace(),  # or whatever is appropriate
  algorithm = "meanfield",
  iter = 40000
)

posterior_samples_fit <- posterior_samples(fit)

h_hat_samples <- posterior_samples_fit$b_h_hat
earning_samples <- posterior_samples_fit$b_earning

diffs <- h_hat_samples - earning_samples

proportion_greater <- mean(diffs > 0)
proportion_smaller <- mean(diffs < 0)

# If proportion_greater is close to 1, then you have strong evidence that 
# h_hat is more strongly associated with happiness than earning is.

n_bootstraps <- 5000
bootstrap_results <- numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  bootstrap_sample <- sample(diffs, size = length(diffs), replace = TRUE)
  bootstrap_results[i] <- mean(bootstrap_sample > 0)
}

lower_bound <- quantile(bootstrap_results, 0.025)
upper_bound <- quantile(bootstrap_results, 0.975)

c(lower_bound, upper_bound)


mod_happiness <- brm(
  h ~ h_hat + reversal * (trial + ema_session) +
    (happiness_hat + reversal * (trial + ema_session) | user_id),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = final_df
)

pp_check(mod_happiness) + xlim(-10, 10)
bayes_R2(mod_happiness)
summary(mod_happiness)
marginal_effects(mod_happiness, "h_hat")

predicted_values <- posterior_predict(mod_happiness, newdata = final_df)

# The output will be a matrix where each row corresponds to an observation in final_df
# and each column corresponds to a posterior draw. You may want to summarize this 
# into a single predicted value per observation, e.g., by taking the mean across columns.

# Take the mean across the columns to get a single predicted value per observation
mean_predicted_values <- colMeans(predicted_values)


final_df |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness)
  )

foo <- final_df |> 
  group_by(reversal, trial) |> 
  summarize(
    h = mean(happiness, trim = 0.1),
    h_hat = mean(happiness_hat, trim = 0.1)
  )
foo$reversal <- as.factor(foo$reversal)

foo |> 
  ggplot(aes(x=trial, y=h, color = reversal)) +
  geom_line()


foo <- final_df |> 
  group_by(reversal, ema_session) |> 
  summarize(
    h = mean(happiness, trim = 0.1),
    h_hat = mean(happiness_hat, trim = 0.1)
  )
foo$reversal <- as.factor(foo$reversal)

foo |> 
  ggplot(aes(x=ema_session, y=h, color = reversal)) +
  geom_line()


