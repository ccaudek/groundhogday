# Analysis of momentary well-being in the PRL task.

# Load necessary library
library(dplyr)



# ------------------------------------------------------------------------------
# Groundhog data

# Function to compute the predicted happiness given parameters and data
compute_predicted_happiness <- function(params, data) {
  w0 <- params[1]
  w1 <- params[2]
  w2 <- params[3]
  w3 <- params[4]
  w4 <- params[5]
  gamma <- params[6]
  
  predicted_happiness <- numeric(nrow(data))
  
  for (t in 1:nrow(data)) {
    predicted_happiness[t] <- w0 +
      w1 * sum(gamma^(t-1:(t-1)) * data$outcome[1:(t-1)]) +
      w2 * sum(gamma^(t-1:(t-1)) * data$reversal[1:(t-1)]) +
      w3 * sum(gamma^(t-1:(t-1)) * data$stimulus[1:(t-1)]) +
      w4 * sum(gamma^(t-1:(t-1)) * data$delta_p[1:(t-1)])
  }
  
  return(predicted_happiness)
}

# Negative log-likelihood assuming normally distributed residuals
nll <- function(params, data) {
  predicted_happiness <- compute_predicted_happiness(params, data)
  ssr <- sum((data$happiness - predicted_happiness)^2)
  n <- nrow(data)
  nll_value <- n/2 * log(2 * pi) + n/2 * log(ssr/n) + n/2
  return(nll_value)
}

# Set seed for reproducibility
set.seed(123)


# Standardize instant mood by user_id.
dz <- d %>%
  group_by(user_id) %>%
  mutate(
    zim = as.vector(scale(instant_mood))
    ) %>%
  ungroup()

# Remove participants because of convergence problems.
bad_id_indices <- c(93, 99, 127, 146, 195, 216)
user_id_codes <- unique(dz$user_id)

bad_codes <- c(
  user_id_codes[bad_id_indices], 
  3338029881, 3665345709, 3248648540
)

good_codes <- setdiff(user_id_codes, bad_codes)
dz_clean <- dz[dz$user_id %in% good_codes, ]
length(unique(dz_clean$user_id))



# Define function to process a single user
process_user <- function(id) {
  
  set.seed(12345)
  
  onesubj_data <- dz_clean |> 
    dplyr::filter(user_id == id)
  
  n_ema_episodes <- length(unique(onesubj_data$ema_number))
  
  par_list <- list()
  
  for (i in seq_len(n_ema_episodes)) {
    
    ema_session <- onesubj_data |> 
      dplyr::filter(ema_number == i)
    
    df <- data.frame(
      trial = ema_session$trial,
      stimulus = ifelse(
        ema_session$is_target_chosen == 0, -1, ema_session$is_target_chosen
      ), # 1 for stimulus A, -1 for stimulus B
      reversal = c(rep(0, 14), 1, rep(0, 15)),  # Reversal occurs at trial 15
      outcome = ifelse(ema_session$feedback == 0, -1, ema_session$feedback),
      delta_p = c(rep(0, 14), 0.6, rep(0, 15)),  # Change in probability at reversal
      happiness = ema_session$zim # standardized by user_id
    )
    
    # Optimize
    init_params <- c(0, 0, 0, 0, 0, 0.5)  # Initial guesses for w0, w1, w2, w3, w4, and gamma
    opt_result <- optim(init_params, nll, data=df)
    mle_params <- opt_result$par
    # add further information
    out <- c(
      mle_params, 
      ifelse(unique(ema_session$is_reversal) == "yes", 1, 0), i, id
    )
    
    par_list[[i]] <- out
  }
  
  # Convert list to dataframe
  par_df <- do.call(rbind, par_list)
  par_df <- as.data.frame(par_df)
  colnames(par_df) <- c(
    "w0", "w1", "w2", "w3", "w4", "gamma", 
    "is_reversal", "ema_number", "user_id"
  )
  
  cat('user_id:', unique(onesubj_data$user_id), '\n')
  
  return(par_df)
}

# Get list of unique user_ids
user_id_codes <- unique(dz_clean$user_id)
length(user_id_codes)

# Apply process_user function to each user_id
results_list <- lapply(user_id_codes, process_user)

# Bind all data frames together into a single data frame
all_results_df <- bind_rows(results_list)

bysubj_mood_df <- dz_clean |> 
  group_by(user_id, ema_number) |> 
  summarize(
    mood_pre = mean(mood_pre),
    mood_post = mean(mood_post),
    control = mean(control)
  ) |> 
  ungroup()


results_df <- left_join(
  all_results_df, bysubj_mood_df, by = c("user_id", "ema_number")
  )



# brms -------------------------------------------------------------------------

plot(density(results_df$w0))
plot(density(results_df$w1))
plot(density(results_df$w2))
plot(density(results_df$w3))
plot(density(results_df$w4))
plot(density(results_df$gamma))

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


results_df$cw1 <- remove_outliers(results_df$w1)
results_df$cw2 <- remove_outliers(results_df$w2)
results_df$cw3 <- remove_outliers(results_df$w3)
results_df$cw4 <- remove_outliers(results_df$w4)
results_df$cwg <- remove_outliers(results_df$gamma)
results_df$cmood_pre <- ifelse(
  results_df$mood_pre == 49 | results_df$mood_pre == -49, NA, 
  results_df$mood_pre)
results_df$cmood_post <- ifelse(
  results_df$mood_post == 49 | results_df$mood_post == -49, NA, 
  results_df$mood_post)

results_df$cmood_pre <- remove_outliers(results_df$cmood_pre)
results_df$cmood_post <- remove_outliers(results_df$cmood_post)


results_w_df <- results_df |> 
  dplyr::select(!c(w1, w2, w3, w4, gamma, mood_pre, mood_post))

imputed_cart = complete(mice(results_w_df, method = "cart"))

# Outcome
imputed_cart$is_reversal <- factor(imputed_cart$is_reversal)
imputed_cart$rev <- ifelse(imputed_cart$is_reversal == 1, 1, -1)

imputed_cart$ema_number_c <- 
  imputed_cart$ema_number - mean(imputed_cart$ema_number)

# Center mood_post by user_id.
bysubj_df <- imputed_cart |> 
  group_by(user_id) |> 
  mutate(
    mood_post_z = as.vector(scale(cmood_post, center = TRUE, scale = TRUE)),
    mood_pre_z = as.vector(scale(cmood_pre, center = TRUE, scale = TRUE)),
    control_z = as.vector(scale(control, center = TRUE, scale = TRUE)),
    # Approximately centered: the fifth session is considered the reference
    ema_number_c = ema_number - 5,
    # difference on the row data. The outliers on mood_pre and mood_post have 
    # been replaced with NAs and then imputed.
    mood_dif = cmood_post - cmood_pre
  ) |> 
  ungroup()




mod <- brm(
  mood_dif ~ cmood_pre + 
    is_reversal * (cw1 + cw2 + cw3 + cw4 + cwg) +
    (is_reversal | user_id/ema_number_c),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = bysubj_df
)

pp_check(mod) 
bayes_R2(mod)
summary(mod)
conditional_effects(mod, "cmood_pre:is_reversal")











mod <- brm(
  dif ~ pre + is_reversal * (w1 + w2 + w3 + w4 + gamma) + 
    (1 + is_reversal | user_id),
  family = asym_laplace(),
  algorithm = "meanfield",
  iter = 40000,
  data = dd
)


bysubj_long_df <- imputed_cart %>%
  pivot_longer(
    cols = starts_with("cmood_"),
    names_to = "cond",
    values_to = "mood",
    values_drop_na = TRUE
  )

fm <- brm(
  mood ~ cond * is_reversal +
    (cond * is_reversal | user_id/ema_number),
  family = asym_laplace(),
  algorithm = "meanfield",
  iter = 40000,
  data = bysubj_long_df
)pp_check(mod) 
summary(mod)
bayes_R2(mod)
conditional_effects(mod, "z1:is_reversal")



hist(results_df$mood_post)

temp <- results_df
temp$z0 <- as.vector(scale(temp$w0))
temp$z1 <- as.vector(scale(temp$w1))
temp$z2 <- as.vector(scale(temp$w2))
temp$z3 <- as.vector(scale(temp$w3))
temp$z4 <- as.vector(scale(temp$w4))
temp$zg <- as.vector(scale(temp$gamma))
temp$mood_post <- as.vector(scale(temp$mood_post))
temp$control <- as.vector(scale(temp$control))
temp$mood_pre <- as.vector(scale(temp$mood_pre))

# Center mood_post by user_id.
bysubj_df <- results_df |> 
  group_by(user_id) |> 
  mutate(
    mood_post_c = as.vector(scale(mood_post, center = TRUE, scale = FALSE)),
    mood_pre_z = as.vector(scale(mood_pre, center = TRUE, scale = TRUE)),
    control_z = as.vector(scale(control, center = TRUE, scale = TRUE)),
    mood_dif = mood_post - mood_pre
  ) |> 
  ungroup()

only_rev_df <- temp |> 
  dplyr::filter(is_reversal != "1")


mod <- brm(
  mood_dif ~ is_reversal * (w1 + w2 + w3 + w4 + gamma) + 
    control_z + mood_pre_z +
    (is_reversal | user_id/ema_number),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = bysubj_df
)
pp_check(mod) + xlim(-4, 4)
summary(mod)
bayes_R2(mod)
conditional_effects(mod, "z1:is_reversal")
conditional_effects(mod, "z2")
conditional_effects(mod, "z3:is_reversal")
conditional_effects(mod, "z4:is_reversal")
conditional_effects(mod, "zg:is_reversal")
conditional_effects(mod, "control_c")
conditional_effects(mod, "mood_pre_c")



fm1 <- brm(
  w1 ~ is_reversal * ema_number_c + 
    (mood_pre + control + rev * ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
  )
pp_check(fm1) + xlim(-4, 4)
summary(fm1)
conditional_effects(fm1, "ema_number_c:rev")

# Reversal
fm2 <- brm(
  w2 ~ rev * ema_number_c + 
    (mood_pre + control + rev * ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm2) + xlim(-4, 4)
summary(fm2)
marginal_effects(fm2, "mood_pre")


# Stimulus
fm3 <- brm(
  w3 ~ rev * ema_number_c + (rev * ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm3) + xlim(-4, 4)
summary(fm3)
marginal_effects(fm3, "ema_number_c:rev")

# Delta P
fm4 <- brm(
  w4 ~ rev * ema_number_c + (rev * ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm4) + xlim(-8, 4)
summary(fm4)
marginal_effects(fm4, "ema_number_c:rev")


# gamma
fm5 <- brm(
  gamma ~ rev * ema_number_c + (rev * ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm5) + xlim(-1, 1)
summary(fm5)
conditional_effects(fm5, "ema_number_c:rev")




# eof ----


all_results_df <- all_results_df |> 
  dplyr::filter(is_reversal == 1)


hist(all_results_df$w1)
hist(all_results_df$w2)
hist(all_results_df$w3)
hist(all_results_df$w4)
hist(all_results_df$gamma)

# Outcome
all_results_df$ema_number_c <- 
  all_results_df$ema_number - mean(all_results_df$ema_number)


fm1 <- brm(
  w1 ~ ema_number_c + (ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm1) + xlim(-4, 4)
summary(fm1)
conditional_effects(fm1, "ema_number_c")

# Reversal
fm2 <- brm(
  w2 ~ ema_number_c + (ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm2) + xlim(-4, 4)
summary(fm2)
conditional_effects(fm2, "ema_number_c")


# Stimulus
fm3 <- brm(
  w3 ~ ema_number_c + (ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm3) + xlim(-4, 4)
summary(fm3)
conditional_effects(fm3, "ema_number_c")

# Delta P
fm4 <- brm(
  w4 ~ ema_number_c + (ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm4) + xlim(-8, 4)
summary(fm4)
conditional_effects(fm4, "ema_number_c")


# gamma
fm5 <- brm(
  gamma ~ ema_number_c + (ema_number_c | user_id),
  family = student(),
  data = all_results_df,
  algorithm = "meanfield"
)
pp_check(fm5) + xlim(-1, 1)
summary(fm5)
conditional_effects(fm5, "ema_number_c")






