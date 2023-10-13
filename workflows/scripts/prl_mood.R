# Script name: prl_mood.R
# Project: groundhog_day
# Script purpose: association between accuracy and mood
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Jun  6 11:13:37 2023
# Last Modified Date: Fri Jun  9 08:52:17 2023
#
# 👉 
#   https://www.frontiersin.org/articles/10.3389/fpsyg.2019.01698/full


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("sjPlot")
library("sjstats")
library("lme4")
library("brms")
library(bayesplot)
library("effectsize")
library("scales")
library("sjPlot")
library("sjmisc")



# IMPORT DATA ------------------------------------------------------------------

# Only reversal data
# d1 <- readRDS("data/prep/groundhog_clean.RDS")

d1 <- readRDS("data/prep/groundhog_all_clean.RDS")


# DATA WRANGLING ---------------------------------------------------------------

d1$user_id <- factor(d1$user_id)

d <- d1 |> 
  dplyr::filter(!is.na(user_id) & ema_number < 12)

length(unique(d$user_id))
# [1] 224


# The variable ema_number represents the consecutive order of the EMA 
# sessions for each subject.
unique_ema_by_user <- lapply(split(d$ema_number, d$user_id), unique)
unique_ema_by_user

# Compliance
temp <- d |>
  group_by(user_id) |>
  summarize(
    mx = max(ema_number)
  )

mean(temp$mx)
# [1] 7.647343


# QUESTION ---------------------------------------------------------------------
# Does is_target_chosen depend on mood_pre?

tgt_data <- d |> 
  # dplyr::filter(is_reversal == "yes") |> 
  dplyr::select(
    is_reversal, is_target_chosen, feedback,
    accuracy, mood_pre, control, trial, ema_number, user_id
  )
tgt_data$mood_pre_z <- as.vector(scale(tgt_data$mood_pre))
tgt_data$control_z <- as.vector(scale(tgt_data$control))
tgt_data$trial_c <- 
  as.vector(scale(tgt_data$trial, center = TRUE, scale = FALSE))
tgt_data$ema_number_c <- 
  as.vector(scale(tgt_data$ema_number, center = TRUE, scale = FALSE))


# glmer analysis ---------------------------------------------------------------

mod_l1 <- glmer(
  is_target_chosen ~ mood_pre_z + control_z + trial_c + ema_number_c + 
    (mood_pre_z + control_z + trial_c + ema_number_c | user_id),
  family = binomial(),
  data = tgt_data
)
plot_model(mod_l1, type = "pred", terms = c("trial_c"))
plot_model(mod_l1, type = "pred", terms = c("ema_number_c"))
plot_model(mod_l1, type = "pred", terms = c("control_z"))
plot_model(mod_l1, type = "pred", terms = c("mood_pre_z"))


# brms analysis ----------------------------------------------------------------

mod1 <- brm(
  is_target_chosen ~ is_reversal*mood_pre_z + control_z + trial_c + 
    ema_number_c + 
    (is_reversal*mood_pre_z + control_z + trial_c + ema_number_c | user_id),
  family = bernoulli(),
  algorithm = "meanfield",
  init = 0.1, 
  data = tgt_data
)

pp_check(mod1)

pp_check(mod1,
         type ="stat",
         stat = "mean",
         ndraws = 1000,
         binwidth = .001)

conditional_effects(mod1, "trial_c")
conditional_effects(mod1, "ema_number_c")
conditional_effects(mod1, "control_z")
conditional_effects(mod1, "is_reversal")
conditional_effects(mod1, "mood_pre_z:is_reversal")

bayes_R2(mod1)


# LOO comparison ---------------------------------------------------------------

md1 <- brm(
  is_target_chosen ~ mood_pre_z + (mood_pre_z | user_id),
  family = bernoulli(),
  algorithm = "meanfield",
  init = 0.1, 
  data = tgt_data
)

md1a <- brm(
  is_target_chosen ~ 1 + (1 | user_id),
  family = bernoulli(),
  algorithm = "meanfield",
  init = 0.1, 
  # backend = "cmdstanr",
  data = tgt_data
)

loo1 <- add_criterion(md1, "loo")
loo1a <- add_criterion(md1a, "loo")
loo_compare(loo1, loo1a)
#       elpd_diff se_diff
# loo1      0.0       0.0
# loo1a -1692.1      58.4

# INTERPRETATION:
# The choice of the target has the highest probability of reward in 
# the first 15 trials, but the lowest probability of reward in the last 15 
# trials. This choice is affected by the present mood at the beginning of the
# session: target choice increases for more positive mood.


# QUESTION ---------------------------------------------------------------------
# Does feedback depend on mood_pre?


mod_l2 <- glmer(
  feedback ~ mood_pre_z + control_z + trial_c + ema_number_c + 
    (mood_pre_z + control_z + trial_c + ema_number_c | user_id),
  family = binomial(link = "logit"),
  data = tgt_data
)
plot_model(mod_l2, type = "pred", terms = c("trial_c"))
plot_model(mod_l2, type = "pred", terms = c("ema_number_c"))
plot_model(mod_l2, type = "pred", terms = c("control_z"))
plot_model(mod_l2, type = "pred", terms = c("mood_pre_z"))

summary(mod_l2)

norev_df <- tgt_data |> 
  dplyr::filter(is_reversal == "no")

norev_df <- norev_df %>%
  group_by(user_id) %>%
  mutate(days = dense_rank(ema_number))

foo <- norev_df |> 
  group_by(user_id) |> 
  summarize(
    n = n_distinct(days),
  ) 
hist(foo$n)

norev_df$days_c <- norev_df$days - 3

mod2 <- brm(
  feedback ~ mood_pre_z + control_z + trial_c + days_c + 
    (mood_pre_z + control_z + trial_c + days_c | user_id),
  family = bernoulli(link = "logit"),
  algorithm = "meanfield",
  init = 0.01, 
  data = norev_df
)

pp_check(mod2)
summary(mod2)

conditional_effects(mod2, "mood_pre_z")
conditional_effects(mod2, "control_z")
conditional_effects(mod2, "trial_c")
conditional_effects(mod2, "days_c")

rev_df <- tgt_data |> 
  dplyr::filter(is_reversal == "yes")

rev_df <- rev_df %>%
  group_by(user_id) %>%
  mutate(days = dense_rank(ema_number))

foo <- rev_df |> 
  group_by(user_id) |> 
  summarize(
    n = n_distinct(days),
  ) 
hist(foo$n)

rev_df <- rev_df |> 
  dplyr::filter(days < 11)

rev_df$days_c <- rev_df$days - 5

mod3 <- brm(
  feedback ~ mood_pre_z + control_z + trial_c + days_c + 
    (mood_pre_z + control_z + trial_c + days_c | user_id),
  family = bernoulli(link = "logit"),
  algorithm = "meanfield",
  init = 0.01, 
  data = norev_df
)

pp_check(mod2)
summary(mod2)

conditional_effects(mod3, "mood_pre_z")
conditional_effects(mod3, "control_z")
conditional_effects(mod3, "trial_c")
conditional_effects(mod3, "days_c")

feedback_df <- bind_rows(rev_df, norev_df)

mod4 <- brm(
  feedback ~ is_reversal* (mood_pre_z + control_z + trial_c + days_c) + 
    (is_reversal*(mood_pre_z + control_z + trial_c + days_c) | user_id),
  family = bernoulli(link = "logit"),
  algorithm = "meanfield",
  init = 0.01, 
  data = feedback_df
)

conditional_effects(mod4, "mood_pre_z:is_reversal")
conditional_effects(mod4, "control_z:is_reversal")
conditional_effects(mod4, "trial_c:is_reversal")
conditional_effects(mod4, "days_c:is_reversal")

bysubj_feedback_df <- feedback_df |> 
  group_by(user_id, is_reversal, days) |> 
  summarize(
    mood_pre_z = mean(mood_pre_z),
    pos_feedback = mean(feedback)
  )

mod5 <- brm(
  pos_feedback ~ is_reversal * mood_pre_z + days + 
    (is_reversal * mood_pre_z + days | user_id),
  family = Beta(),
  algorithm = "meanfield",
  init = 0.01, 
  data = bysubj_feedback_df
)

conditional_effects(mod5, "mood_pre_z:is_reversal")
bayes_R2(mod5)

# INTERPRETATION:
# The results on feedback are similar to those on is_target_chosen.


# QUESTION ---------------------------------------------------------------------
# Does mood change (post - pre) depend on ema_number?

mood_df <- d |> 
  # dplyr::filter(mood_post != 0 & mood_pre != 0) |> 
  group_by(user_id, ema_number) |> 
  summarize(
    mood_change = mean(mood_change, na.rm = TRUE),
    mood_post = mean(mood_post, na.rm = TRUE),
    mood_pre = mean(mood_pre, na.rm = TRUE),
    control = mean(control, na.rm = TRUE),
    gain = mean(gain, na.rm = TRUE)
  ) |> 
  ungroup()

# mood_change = mood_post - mood_pre
mood_df |> 
  group_by(ema_number) |> 
  summarize(
    mood_change = mean(mood_change, na.rm = TRUE)
  ) 

mood_df$mood_ch_z <- as.vector(scale(mood_df$mood_change)) 
mood_df$mood_pre_z <- as.vector(scale(mood_df$mood_pre)) 
mood_df$mood_post_z <- as.vector(scale(mood_df$mood_post)) 
mood_df$control_z <- as.vector(scale(mood_df$control)) 
mood_df$ema_number_z <- as.vector(scale(mood_df$ema_number)) 
mood_df$gain_z <- as.vector(scale(mood_df$gain)) 


mod_l3 <- lmer(
  mood_ch_z ~ gain_z + ema_number_z + mood_pre_z + control_z +
    (ema_number_z + mood_pre_z + control_z | user_id),
  data = mood_df
)
summary(mod_l3)


mod3 <- brm(
  mood_ch_z ~ gain_z + ema_number_z + mood_pre_z + 
    (gain_z + ema_number_z + mood_pre_z | user_id),
  family = asym_laplace(),
  algorithm = "meanfield",
  data = mood_df
)

# simulations vs. obs: Overall
pp_check(mod3,
         type ="dens_overlay",
         ndraws = 100)

pp_check(mod3,
         type ="stat",
         stat = "mean",
         ndraws = 1000,
         binwidth = .001)


pp_check(mod3,
         type ="stat",
         stat = "sd",
         ndraws = 500,
         binwidth = .001)


pp_check(mod3,
         type ="stat_2d")


pp_check(mod3,
         type ="scatter_avg")



# MODEL EVAL: LOO --------------------------------------------------------

# loo & pareto K
model_loo <- loo(mod3, save_psis = TRUE, cores = 4)
plot(model_loo, diagnostic = "k")
plot(model_loo, diagnostic = "n_eff")

# loo pit
w <- weights(model_loo$psis_object)
ppc_loo_pit_overlay(y = mood_df$mood_ch_z, 
                    yrep = posterior_predict(mod3), 
                    lw = w)
ppc_loo_pit_qq(y = mood_df$mood_ch_z, 
               yrep = posterior_predict(mod3), 
               lw = w)

# MODEL EVAL: COND. EFF. -------------------------------------------------------

plot(
  conditional_effects(
    mod3, 
    "gain_z",
    re_formula = NA,
    method = "fitted",
    spaghetti = TRUE,
    ndraws = 200
  ),
  points = FALSE
)

plot(
  conditional_effects(
    mod3, 
    "ema_number_z",
    re_formula = NA,
    method = "fitted",
    spaghetti = TRUE,
    ndraws = 200
  ),
  points = FALSE
)


plot(
  conditional_effects(
    mod3, 
    "mood_pre_z",
    re_formula = NA,
    method = "fitted",
    spaghetti = TRUE,
    ndraws = 200
  ),
  points = FALSE
)




# QUESTION  
# Does gain increase as a function of ema_number?

d |> 
  group_by(ema_number) %>%
  summarise(
    sdterr = sd(gain, na.rm = TRUE) / sqrt(n()),
    gain = mean(gain)
  ) |> 
  ggplot(aes(x=ema_number, y=gain)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=gain-sdterr, ymax=gain+sdterr), width=.2)


# unique_subjects <- unique(d$user_id)
# num_subjects = 1
# selected_subjects <- sample(unique_subjects, size = num_subjects, replace = FALSE)
# selected_data <- d[d$user_id %in% selected_subjects, ]  
# selected_data$epoch <- factor(selected_data$epoch)
# 
# df_for_plot <- selected_data |> 
#   group_by(user_id, epoch, ema_number) |> 
#   summarize(
#     mood_post = mean(mood_post, na.rm = TRUE),
#     gain = mean(gain, na.rm = TRUE),
#     TIME_total = mean(TIME_total, na.rm = TRUE)
#   ) |> 
#   ungroup()
# 
# # Lattice plot for mood_post vs. days
# df_for_plot |> 
#   ggplot(aes(x=factor(TIME_total),y=gain, group=1)) + 
#   # geom_line() +
#   geom_point() +
#   geom_smooth(method=lm, se=FALSE) +
#   facet_wrap(~user_id, ncol=5) +   
#   theme(strip.text.x=element_blank()) + 
#   labs(x="EMA number",y="Affect")

# 
# df_bysubj <- d |> 
#   group_by(user_id, ema_number) |> 
#   summarize(
#     mood_pre = mean(mood_pre, na.rm = TRUE),
#     mood_post = mean(mood_post, na.rm = TRUE),
#     gain = mean(gain, na.rm = TRUE),
#     accuracy = mean(accuracy, na.rm = TRUE),
#     TIME_total = mean(TIME_total, trim = 0.1,  na.rm = TRUE)
#   ) |> 
#   mutate(mood_post = rescale(mood_post, to = c(-1, 1))) %>%
#   ungroup()
# 
# 
# fm <- lmer(
#   gain ~ mood_post + ema_number +
#     (1 + mood_post + ema_number | user_id), 
#   control = lmerControl(optimizer ="Nelder_Mead"),
#   # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
#   df_bysubj
# )
# summary(fm)
# plot_model(fm, "eff", terms = c("mood_post"))
# plot_model(fm, "eff", terms = c("ema_number"))

# 
# fm <- lmer(
#   mood_post ~ gain + ema_number + 
#     (1 + gain + ema_number | user_id), 
#   # control = lmerControl(optimizer ="Nelder_Mead"),
#   control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
#   df_bysubj
# )
# summary(fm)
# summary(rePCA(fm))
# 
# plot_model(fm, "eff", terms = c("days", "epoch", "mood_post"))
# 
# m1 <- brm(
#   gain ~ mood_post * ema_number + 
#     (1 + mood_post + ema_number | user_id), 
#   family = asym_laplace(),
#   backend = "cmdstanr",
#   data = df_bysubj
# )
# 
# summary(m1)
# conditional_effects(m1, "ema_number:mood_post")
# pp_check(m1)
# 
# 
# hist((df_bysubj$gain))



# QUESTION
# Instant mood before and after the threshold


d |> 
  group_by(user_id, epoch, trial) %>%
  summarise(
    im = mean(instant_mood),
  ) |> 
  group_by(trial) %>%
  summarise(
    imood = mean(im),
    se = sqrt(var(im, na.rm = TRUE) / n())
  ) |> 
  ggplot(aes(x=trial, y=imood)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=imood-se, ymax=imood+se), width=.2)




bysubj_df <- d |> 
  group_by(user_id, epoch, trial) %>%
  summarise(
    im = mean(instant_mood),
  )

bysubj_df$im |> hist()
bysubj_df$epoch <- factor(bysubj_df$epoch)

mod_1 <- brm(
  im ~ trial * epoch + (1 + trial * epoch | user_id), 
  family = student(),
  # backend="cmdstanr",
  algorithm = "meanfield",
  data = bysubj_df
)
summary(mod_1)
pp_check(mod_1)
conditional_effects(mod_1, "trial:epoch")


m2 <- brm(
  gain ~ (mood_pre + mood_post) * days +
    (1 + mood_pre + mood_post + days | user_id),
  family = gaussian(),
  backend = "cmdstanr",
  data = df_by_subj
)

summary(m2)
conditional_effects(m2, "mood_post")
conditional_effects(m2, "days")
performance::r2_bayes(m2)
performance::r2_loo(m2)


hist(d$instant_mood)

foo <- d |> 
  group_by(user_id, ema_number) |> 
  summarize(
    instant_mood = mean(instant_mood),
    mood_pre = mean(mood_pre),
    mood_post = mean(mood_post),
    gain = mean(gain),
    TIME_total = mean(TIME_total)
  ) |> 
  ungroup()

fm <- lmer(
  instant_mood ~ 1 + ema_number + mood_pre +
    (1 + ema_number | user_id),
  foo
)
summary(fm)

# Lattice plot for NA vs. Performance Type
d |> 
  # dplyr::filter(ema_number < 9) |> 
  group_by(trial) |> 
  summarize(
    instant_mood = mean(instant_mood)
  ) |> 
  #dplyr::filter(ema_number == 1) |> 
ggplot(aes(x=factor(trial), y=instant_mood)) + 
  geom_point() + 
  #facet_wrap(~ema_number,ncol=4) +   
  #theme(strip.text.x=element_blank()) +
  labs(x="Instant mood",y="Gain")

foo <- d |> 
  group_by(trial) |> 
  summarize(
    imood = mean(instant_mood)
  )

fm <- lm(
  imood ~ poly(trial, 3),
  foo
)

plot(foo$trial, foo$imood)
lines(sort(foo$trial),                 # Draw polynomial regression curve
      fitted(fm)[order(foo$trial)],
      col = "red",
      type = "l")

#----------------
# Instant mood.
#----------------
# Regression Discontinuity
# https://mixtape.scunning.com/06-regression_discontinuity

m0 <- lmer(
  instant_mood ~ 1 + (1 | user_id) + (1 | ema_number),
  data = d
)
summary(m0)  

icc = 0.25737 / (0.25737 + 0.07843 + 1.05783)
icc
# 18.5% of the total variability in performance anxiety scores are 
# attributable to differences among subjects. In this particular model, 
# we can also say that the average correlation for any pair of responses 
# from the same individual is a moderately low .185.

d$user_day <- paste(d$user_id, d$ema_number, sep="_")

library(splines)
m1 <- lmer(
  instant_mood ~ 1 + ns(trial, 3) + 
    (1 + ns(trial, 3) | user_id) + (1 | user_day),
  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
  data = d
)
summary(m1) 

d$predicted <- predict(m1, newdata = d, type = "response")

df_forplot <- d |> 
  group_by(trial) |> 
  summarize(
    y = mean(instant_mood),
    yhat = mean(predicted)
  )

plot(df_forplot$trial, df_forplot$y, pch = 16, xlab = "Predictor 1", ylab = "Judgment")
lines(df_forplot$trial, df_forplot$yhat, col = "blue", lwd = 2)



m3 <- lmer(
  instant_mood ~ 1 + trial + 
    (1 + trial | user_id) + (1 | user_day),
  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
  data = d, subset=trial < 16
)
summary(m3) 

m4 <- lmer(
  instant_mood ~ 1 + trial + 
    (1 + trial | user_id) + (1 | user_day),
  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
  data = d, subset=trial > 15
)
summary(m4) 


df_byday <- d |> 
  group_by(user_id, ema_number) |> 
  summarize(
    control = mean(control),
    mood_pre = mean(mood_pre),
    mood_post = mean(mood_post)
  ) |> 
  ungroup()
df_byday$user_day <- paste(df_byday$user_id, df_byday$ema_number, sep="_")

bm5 <- brm(
  control ~ 1 + mood_pre + mood_post +
    (1 + mood_pre + mood_post | user_id) + (1 | ema_number),
  family = cumulative("probit"),
  data = df_byday,
  init = 0.1,
  backend = "cmdstanr"
)
pp_check(bm5)
conditional_effects(bm5, "mood_pre", categorical = TRUE)
conditional_effects(bm5, "mood_post", categorical = TRUE)
summary(bm5)

plot(density(df_byday$mood_pre))

bm6 <- brm(
  mood_pre ~ 1 + control +
    (1 + control | user_id) + (1 | ema_number),
  family = skew_normal(),
  data = df_byday,
  init = 0.1,
  backend = "cmdstanr"
)
pp_check(bm6)
conditional_effects(bm6, "control")
summary(bm6)



a <- summary(bm6)
summary_mod1 <- rbind(data.frame(a$fixed), data.frame(a$spec_pars) )
rownames(summary_mod1) <- c("$\\beta_0$", "$\\beta$", "$\\sigma_{e}$", "$\\alpha$")
colnames(summary_mod1) <- c("mean","SE", "lower bound", "upper bound", "Rhat", "Bulk ESS", "Tail ESS")

summary_mod1 %<>%
  select(-c("Bulk ESS", "Tail ESS")) %>% # removing ESS
  rownames_to_column(var = "parameter") 

summary_mod1 <- summary_mod1 |> mutate_if(is.numeric, round, 3)

write.csv(summary_mod1, file = snakemake@output[["table_1"]], row.names = FALSE)






m5 <- lmer(
  control ~ 1 + mood_pre + mood_post +
    (1 + mood_pre + mood_post | user_id) + (1 | ema_number),
  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
  data = df_byday
)
summary(m5) 
summary(rePCA(m5))

fm1 <- lmer(
  instant_mood ~ poly(trial, 3) +
    (poly(trial, 3) | user_id) +
    (1 | ema_number),
  data = d
)
summary(fm1)  


foo <- d |> 
  group_by(user_id, ema_number) |> 
  summarize(
    imood = mean(instant_mood),
    control = mean(control),
    mood_post = mean(mood_post),
    mood_pre = mean(mood_pre),
    gain = mean(gain),
    time = mean(TIME_total)
  ) |> 
  ungroup()

fm2 <- lmer(
  control ~ 1 + mood_pre +
    (1 + mood_pre | user_id) +  (1 | ema_number),
  data = foo
)
summary(fm2)  

library(rdrobust)

foo <- d |> 
  group_by(user_id, trial) |> 
  summarize(
    imood = mean(instant_mood),
  ) |> 
  ungroup()

rdr <- rdrobust(
  y = foo$imood,
  x = foo$trial, 
  c = 16
)
summary(rdr)


fm2 <- lmer(
  gain ~ 1 + control * mood_pre +
    (1 + control + mood_pre | user_id) + (1 | ema_number),
  data = foo
)
summary(fm2)  
