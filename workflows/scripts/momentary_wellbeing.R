# Analysis of momentary well-being in the PRL task.
# This script computes the parameters of the momentary happiness model
# for each user_id and each ema_number.

# Load necessary library
library("here")
library("tidyverse")
library("mice")
library("lme4")
library("brms")
library("bayesplot")
library("effectsize")
library("scales")
library("sjstats")
library("sjPlot")
library("sjmisc")

# Functions:
source(here::here("workflows", "scripts", "funs", "funs_instant_mood.R"))


# This script requires the dataframe d from prl_mood.R.
d1 <- readRDS("data/prep/groundhog_all_clean.RDS")

# DATA WRANGLING
d1$user_id <- as.numeric(d1$user_id)
# Remove last ema sessions with too low compliance
d <- d1 |> 
  dplyr::filter(!is.na(user_id) & ema_number < 13)

length(unique(d$user_id))
# [1] 224


# Data wrangling ---------------------------------------------------------------

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

# Get list of unique user_ids
user_id_codes <- unique(dz_clean$user_id)
length(user_id_codes)

# Apply process_user function to each user_id
par_list <- NULL
results_list <- lapply(user_id_codes, process_user)

# Bind all data frames together into a single data frame
all_results_df <- bind_rows(results_list)

# Add mood_pre, mood_post, control.
bysubj_mood_df <- dz_clean |> 
  group_by(user_id, ema_number) |> 
  summarize(
    mood_pre = mean(mood_pre),
    mood_post = mean(mood_post),
    control = mean(control)
  ) |> 
  ungroup()

# Final dataframe.
results_df <- left_join(
  all_results_df, bysubj_mood_df, by = c("user_id", "ema_number")
  )

names(results_df)
# [1] "w0"          "w1"          "w2"          "w3"          "w4"         
# [6] "w5"          "gamma"       "is_reversal" "ema_number"  "user_id"    
# [11] "alpha"       "mood_pre"    "mood_post"   "control"  


# brms -------------------------------------------------------------------------

plot(density(results_df$w0))
plot(density(results_df$w1))
plot(density(results_df$w2))
plot(density(results_df$w3))
plot(density(results_df$w4))
plot(density(results_df$w5))
plot(density(results_df$gamma))
plot(density(results_df$alpha))


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
results_df$cw5 <- remove_outliers(results_df$w5)
results_df$cwg <- remove_outliers(results_df$gamma)

results_df$cmood_pre <- remove_outliers(results_df$mood_pre)
results_df$cmood_post <- remove_outliers(results_df$mood_post)

FLAG_MOOD <- 0 # 0: keep mood data unchanged; 1: impute the values == 0
if (FLAG_MOOD) {
  # Remove and impute the values == 0
  results_df$cmood_post <- ifelse(
    results_df$cmood_post == 0, NA, results_df$cmood_post
  )

  results_df$cmood_pre <- ifelse(
    results_df$cmood_pre == 0, NA, results_df$cmood_pre
  )
} else {
  results_df$cmood_pre <- results_df$mood_pre
  results_df$cmood_post <- results_df$mood_post
}

# Remove and impute the values alpha = 0 ad alpha = 1
results_df$calpha <- ifelse(
  results_df$alpha == 0 | results_df$alpha == 1, NA, results_df$alpha
)

# Remove redundant columns
results_w_df <- results_df |> 
  dplyr::select(!c(w0, w1, w2, w3, w4, w5, gamma, mood_pre, mood_post))

# Imputation of missing data.
imputed_cart = complete(mice(results_w_df, method = "cart"))

# Data wrangling.
imputed_cart$is_reversal <- factor(imputed_cart$is_reversal)
imputed_cart$rev <- ifelse(imputed_cart$is_reversal == 1, 1, -1)

imputed_cart$ema_number_c <- 
  imputed_cart$ema_number - mean(imputed_cart$ema_number)

# Difference post - pre on the raw data. 
bysubj_df <- imputed_cart |> 
  # group_by(user_id) |> 
  mutate(
    # Approximately centered: the fifth session is considered the reference
    ema_number_c = ema_number - 5,
    # The outliers on mood_pre and mood_post have been replaced with NAs and 
    # then imputed.
    mood_dif = cmood_post - cmood_pre
  ) |> 
  dplyr::rename(
    mood_pre = cmood_pre,
    mood_post = cmood_post  
  ) |> 
  ungroup()
  
bysubj_df$environment <- 
  ifelse(bysubj_df$is_reversal == 1, "Volatile", "Stable") |> 
  as.factor()


# Effect of environment on alpha -----------------------------------------------

bysubj_df$calpha |> hist()

m <- brm(
  calpha ~ environment +
    (environment | user_id / ema_number_c),
  family = asym_laplace(),
  algorithm = "meanfield",
  data = bysubj_df
)
pp_check(m)
summary(m)
bayes_R2(m)
conditional_effects(m, "environment")

# In a highly volatile environment, an agent might benefit from a higher α 
# value, allowing it to rapidly update its expectations to adapt to the 
# changing contingencies. A high α value places more weight on recent outcomes, 
# making the agent more responsive to changes.


# Procedure effect on mood difference ------------------------------------------

bysubj_df$mood_dif |> hist()

mod1 <- brm(
  mood_dif ~ mood_pre + control + 
    environment * (cw1 + cw2 + cw3 + cw4 + cw5 + cwg) +
    (mood_pre + control + environment | user_id / ema_number_c),
  family = student(),
  algorithm = "meanfield",
  data = bysubj_df
)
loo1 <- loo(mod1)
plot(loo1)

pp_check(mod1)  + xlim(-80, 80)
bayes_R2(mod1)
summary(mod1)
# w1: Outcome
# w2: Reversal
conditional_effects(mod1, "cw3:environment") # w3: Stimulus
conditional_effects(mod1, "cw4") # w4: delta_p
conditional_effects(mod1, "cw5:environment") # w5: RPE
conditional_effects(mod1, "cwg:environment") # gamma
conditional_effects(mod1, "environment")
conditional_effects(mod1, "control")
conditional_effects(mod1, "mood_pre")


mod2 <- brm(
  mood_dif ~ mood_pre + control + environment +
    (mood_pre + control + environment | user_id / ema_number_c),
  family = student(),
  algorithm = "meanfield",
  data = bysubj_df
)
loo2 <- loo(mod2)
plot(loo2)

# compare both models
loo_compare(loo1, loo2)  









bysubj_df |> 
  group_by(user_id, environment) |> 
  summarize(
    md = mean(mood_dif, trim = 0.1),
    pre = mean(mood_pre, trim = 0.1),
    post = mean(mood_post, trim = 0.1)
  ) |> 
  group_by(environment) |> 
  summarize(
    mood_dif = mean(md, trim = 0.1),
    sdterr = sqrt(var(md) / n()),
    n = n(),
    mood_pre = mean(pre),
    mood_post = mean(post)
  )

mod2 <- brm(
  mood_dif ~ mood_pre * environment + 
    (mood_pre * environment | user_id / ema_number),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = bysubj_df
)

pp_check(mod2) + xlim(-80, 80)
summary(mod2)
bayes_R2(mod2)
conditional_effects(mod2, "environment")


mod3 <- brm(
  mood_dif ~ control + ema_number + mood_pre + 
    (control + ema_number + mood_pre | user_id),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = bysubj_df
)

pp_check(mod3) 
summary(mod3)
bayes_R2(mod3)
conditional_effects(mod3, "ema_number")
conditional_effects(mod3, "mood_pre")






bysubj_long_df <- bysubj_df |> 
  dplyr::select(-mood_dif) |> 
  pivot_longer(
    cols = starts_with("mood_"),
    names_to = "cond",
    values_to = "mood",
    values_drop_na = TRUE
  )

bysubj_long_df |> 
  group_by(cond) |> 
  summarize(
    avg = mean(mood)
  )

mod3 <- brm(
  mood ~ control + environment * cond +
    (control + environment * cond | user_id/ema_number),
  family = student(),
  algorithm = "meanfield",
  iter = 40000,
  data = bysubj_long_df
)

pp_check(mod3) 
summary(mod3)
bayes_R2(mod3)
conditional_effects(mod3, "environment:cond")



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





