###########################################################################
####                          Breath Duration                          ####
###########################################################################

# Author:         Martin Grund, Marc Pabst
# Last update:    August 6, 2021


# Settings ---------------------------------------------------------------

rm(list=ls())

data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

# Create colormap
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=4
colmap <- c(rgb(141,211,199, maxColorValue=255),
            rgb(255,255,179, maxColorValue=255),
            rgb(190,186,218, maxColorValue=255),
            rgb(251,128,114, maxColorValue=255))


# Load packages
library("cowplot")
library("ggpubr")
library("rstatix")
library("tidyverse")
library(circular)


# Load data ---------------------------------------------------------------

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210622.Rds", sep="/"))


# Filter data -------------------------------------------------------------

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_exhale_filter == 1)
data_stims_valid$resp_cycle_t <- data_stims_valid$resp_cycle_t_exhale
data_stims_valid$stim_degree <- data_stims_valid$stim_degree_exhale

#data_stims_valid$resp_cycle_t <- data_stims_valid$resp_cycle_t_inhale

# General statistics ------------------------------------------------------

# Number of trials with valid respiratory cycle per participants
num_breath <- aggregate(resp_cycle_filter ~ ID, data_stims_valid, FUN = sum)

range(num_breath$resp_cycle_filter)

mean(num_breath$resp_cycle_filter)

# Mean duration per participant
#resp_cycle_t <- aggregate(resp_cycle_t ~ ID, subset(data_stims, resp_cycle_filter == 1), FUN = mean)
resp_cycle_t <- aggregate(resp_cycle_t ~ ID, data_stims_valid, FUN = mean)

mean(resp_cycle_t$resp_cycle_t)
sd(resp_cycle_t$resp_cycle_t)
range(resp_cycle_t$resp_cycle_t)


# Plot CR, miss, and hit --------------------------------------------------

# BY TRIAL TYPE

# Subset Data
resp_cycle_t_data = data_stims_valid %>%
  filter(trial_type != "FA") %>%
  #order
  mutate(trial_type = factor(trial_type, levels=c("FA", "CR", "near_miss", "near_hit"))) %>%
  mutate(trial_type = droplevels(trial_type)) %>%
  group_by(ID, trial_type) %>%
 # dplyr::summarize(mean_resp_cycle_t = mean(resp_cycle_t, na.rm=T)) %>%
  dplyr::summarize(mean_resp_cycle_t = mean(resp_cycle_t, na.rm=T),
                   mean_stim_degree = mean(circular(stim_degree, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0), na.rm=T),
                   var_stim_degree = var(circular(stim_degree, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0), na.rm=T)) %>%
  ungroup()


# model1 <- lmer(mean_resp_cycle_t ~ trial_type + (1|ID), resp_cycle_t_data)
# model2 <- lmer(mean_resp_cycle_t ~ var_stim_degree + (1|ID), resp_cycle_t_data)
# model3 <- lmer(mean_resp_cycle_t ~ trial_type + var_stim_degree + (1|ID), resp_cycle_t_data)


# ANOVA
library("afex")

# # Center mean
# for (c in unique(resp_cycle_t_data$trial_type)) {
#   #print(mean(resp_cycle_t_data$var_stim_degree[resp_cycle_t_data$trial_type == c]))
#   resp_cycle_t_data$var_stim_degree[resp_cycle_t_data$trial_type == c] <- resp_cycle_t_data$var_stim_degree[resp_cycle_t_data$trial_type == c] - mean(resp_cycle_t_data$var_stim_degree[resp_cycle_t_data$trial_type == c])
# }

fit <- aov_ez("ID","mean_resp_cycle_t", resp_cycle_t_data, within=c("trial_type"))
fit # to see corrected degrees of freedom 
summary(fit) # see epsilon values

# Order trial type levels
#revalue(resp_cycle_t_data$trial_type, c('near_miss' = 'Miss', 'near_hit' = 'Hit'))
resp_cycle_t_data$trial_type <- factor(resp_cycle_t_data$trial_type,
                                          levels = c('CR', 'near_miss', 'near_hit'),
                                          ordered = T)

# Output means
aggregate(mean_resp_cycle_t ~ trial_type, resp_cycle_t_data, FUN=mean)
  
# Test Breath Duration
stat_test = resp_cycle_t_data %>% 
  t_test(mean_resp_cycle_t ~ trial_type, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=1.5)

stat_test$groups
stat_test$estimate*1000

# Plot Breath Duration
ggplot(aes(y = mean_resp_cycle_t, x = trial_type), data = resp_cycle_t_data) + 
  geom_boxplot(aes(fill = trial_type)) +
  theme_cowplot() + 
  #ylab("Respiratory cycle duration in s") +
  ggtitle("Mean respiratory cycle duration") +
  stat_pvalue_manual(stat_test, label = "p.adj", tip.length = 0.01, hide.ns=T) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) + 
  scale_fill_manual(values=colmap[c(1,3,4)]) +
  scale_x_discrete(labels = c('CR', 'Miss', 'Hit')) #+
  #scale_y_continuous(breaks = seq(3,5.5,.5))


# ANOVA Detection x Confidence - Near-threshold trials --------------------

data_stims_valid$confidence <- NA
data_stims_valid$confidence[data_stims_valid$resp2==0] <- 'unconfident'
data_stims_valid$confidence[data_stims_valid$resp2==1] <- 'confident'

mean_resp_cycle_t <- aggregate(resp_cycle_t ~ ID*trial_type*confidence, subset(data_stims_valid, trial_type != "FA" & trial_type != "CR"), FUN=mean)

fit <- aov_ez("ID","resp_cycle_t", mean_resp_cycle_t, within=c("trial_type", "confidence"))
fit

mean_resp_cycle_t$cond <- paste(mean_resp_cycle_t$trial_type, mean_resp_cycle_t$confidence, sep = '_')

stat_test = mean_resp_cycle_t %>% 
  t_test(resp_cycle_t ~ cond, paired=TRUE, p.adjust.method = "fdr")

  
# Plot confident/unconfident miss/hit -------------------------------------

# BY CONFIDENCE x TRIAL_TYPE

# Subset Data
resp_cycle_t_data = data_stims_valid %>%
  filter(!is.na( resp_cycle_t ) ) %>%
  filter(trial_type != "FA" & trial_type != "CR") %>%
  #filter(trial_type != "FA") %>%
  mutate(trial_type = paste(trial_type, confidence)) %>%
  group_by(ID, trial_type) %>%
  dplyr::summarize(mean_resp_cycle_t = mean(resp_cycle_t)) %>%
  ungroup()

# resp_cycle_t_data = data_stims_valid %>%
#   filter(!is.na( resp_cycle_t ) ) %>%
#   filter(trial_type != "FA") %>%
# 
# 
# table(mean_resp_cycle_t$ID)
# 
# model1 = lme4::lmer(formula = resp_cycle_t ~ 1 + (1 | ID), data = mean_resp_cycle_t, REML = FALSE)
# model2 = lme4::lmer(formula = resp_cycle_t ~ trial_type + confidence + (1 | ID), data = mean_resp_cycle_t, REML = FALSE)

# ANOVA
library("afex")

fit <- aov_ez("ID","mean_resp_cycle_t", resp_cycle_t_data, within=c("trial_type"))
fit # to see corrected degrees of freedom 
summary(fit) # see epsilon values

# Order trial type levels
resp_cycle_t_data$trial_type <- factor(resp_cycle_t_data$trial_type,
                                          levels = c('near_miss unconfident', 'near_miss confident', 'near_hit unconfident', 'near_hit confident'),
                                          ordered = T)

# Test Breath Duration
stat_test = resp_cycle_t_data %>% 
  t_test(mean_resp_cycle_t ~ trial_type, paired=TRUE, p.adjust.method = "fdr") %>% 
  add_xy_position(step.increase=0.05)

print(stat_test)

# Plot Breath Duration
ggplot(aes(y=mean_resp_cycle_t, x=trial_type), data = resp_cycle_t_data) + 
  geom_boxplot(aes(fill=trial_type)) +
  theme_cowplot() + 
  ylab("Mean respiratory cycle duration") +
  scale_fill_manual(values=colmap[c(3,3,4,4)]) +
  stat_pvalue_manual( stat_test, label = "p.adj", tip.length = 0.01, hide.ns=T) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()) +
  scale_x_discrete(labels = c('Miss C-', 'Miss C+', 'Hit C-', 'Hit C+'))# +
  #scale_y_continuous(breaks = seq(3,6,.5))


# Plot confident and unconfident near-threshold trials ---------------------------------------------

# BY CONFIDENCE

# Subset Data
resp_cycle_t_data = data_stims_valid %>%
  #filter(trial_type != "FA") %>%
  filter(trial_type != "FA" & trial_type != "CR") %>%
  #mutate(trial_type = droplevels(trial_type)) %>%
  mutate(trial_type = trial_type) %>%
  group_by(ID, confidence) %>%
  dplyr::summarize(mean_resp_cycle_t = mean(resp_cycle_t, na.rm=T)) %>%
  ungroup()

# Test Breath Duration
stat_test = resp_cycle_t_data %>% 
  t_test(mean_resp_cycle_t ~ confidence, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

print(stat_test)

# Plot Breath Duration
resp_cycle_t_data %>%
  ggplot(aes(y=mean_resp_cycle_t, x=confidence, color=confidence)) + 
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_brewer(palette="Dark2")  +
  stat_pvalue_manual( stat_test, label = "p", tip.length = 0.01, hide.ns=T)+
  theme(legend.position="none")


# Which phase is adapted --------------------------------------------------

data_test <- subset(data_stims_valid,trial_type != 'FA')

# Exhale duration
mean_exhale_duration_exhale <- aggregate(exhale_duration_exhale ~ ID*trial_type, data_test, FUN = mean)

boxplot(exhale_duration_exhale ~ trial_type, mean_exhale_duration_exhale)

stat_test_exhale = mean_exhale_duration_exhale %>% 
    t_test(exhale_duration_exhale ~ trial_type, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

mean_exhale_duration_exhale %>%
  ggplot(aes(y=exhale_duration_exhale, x=trial_type, color=trial_type)) + 
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_brewer(palette="Dark2")  +
  stat_pvalue_manual( stat_test_exhale, label = "p", tip.length = 0.01, hide.ns=T)+
  theme(legend.position="none")

# Inhale duration
mean_inhale_duration_exhale <- aggregate(inhale_duration_exhale ~ ID*trial_type, data_test, FUN = mean)

boxplot(inhale_duration_exhale ~ trial_type, mean_inhale_duration_exhale)

stat_test_inhale = mean_inhale_duration_exhale %>%
  t_test(inhale_duration_exhale ~ trial_type, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

mean_inhale_duration_exhale %>%
  ggplot(aes(y=inhale_duration_exhale, x=trial_type, color=trial_type)) + 
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_brewer(palette="Dark2")  +
  stat_pvalue_manual( stat_test_inhale, label = "p", tip.length = 0.01, hide.ns=T)+
  theme(legend.position="none")

p.adjust(c(stat_test_exhale$p, stat_test_inhale$p), 'fdr')

# Test respiration cycle duration based on sum
data_test$sum_in_ex_t_exhale <- data_test$exhale_duration_exhale + data_test$inhale_duration_exhale

mean_sum_in_ex_t_exhale <- aggregate(sum_in_ex_t_exhale ~ ID*trial_type, data_test, FUN = mean)

stat_test_sum_in_ex = mean_sum_in_ex_t_exhale %>%
  t_test(sum_in_ex_t_exhale ~ trial_type, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

mean_sum_in_ex_t_exhale %>%
  ggplot(aes(y=sum_in_ex_t_exhale, x=trial_type, color=trial_type)) + 
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_brewer(palette="Dark2")  +
  stat_pvalue_manual( stat_test_sum_in_ex, label = "p", tip.length = 0.01, hide.ns=T)+
  theme(legend.position="none")

# Test respiration cycle duration based on ratio
data_test$ratio_in_ex_t_exhale <- data_test$exhale_duration_exhale/data_test$inhale_duration_exhale

mean_ratio_in_ex_t_exhale <- aggregate(ratio_in_ex_t_exhale ~ ID*trial_type, data_test, FUN = mean)

stat_test_ratio_in_ex = mean_ratio_in_ex_t_exhale %>%
  t_test(ratio_in_ex_t_exhale ~ trial_type, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

mean_ratio_in_ex_t_exhale %>%
  ggplot(aes(y=ratio_in_ex_t_exhale, x=trial_type, color=trial_type)) + 
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_brewer(palette="Dark2")  +
  stat_pvalue_manual( stat_test_ratio_in_ex, label = "p", tip.length = 0.01, hide.ns=T)+
  theme(legend.position="none")


# ANOVA -------------------------------------------------------------------

# Requires mean_exhale_duration_exhale & mean_inhale_duration_exhale from previous section

mean_ex_in_t_exhale <- mean_exhale_duration_exhale

mean_ex_in_t_exhale$inhale_duration_exhale <- NA

for (i in c(1:nrow(mean_ex_in_t_exhale))) {
  mean_ex_in_t_exhale$inhale_duration_exhale[i] <- mean_inhale_duration_exhale$inhale_duration_exhale[mean_inhale_duration_exhale$ID == mean_ex_in_t_exhale$ID[i] & mean_inhale_duration_exhale$trial_type == mean_ex_in_t_exhale$trial_type[i]]
}

colnames(mean_ex_in_t_exhale) <- c('ID', 'trial_type', 'exhale', 'inhale')

mean_ex_in_t <- pivot_longer(mean_ex_in_t_exhale,
                             cols = c(exhale, inhale), 
                             names_to = 'phase', values_to = 'duration')

mean_ex_in_t$trial_type <- as.factor(mean_ex_in_t$trial_type)
mean_ex_in_t$phase <- as.factor(mean_ex_in_t$phase)


library("afex")

fit_all <- aov_ez("ID","duration", mean_ex_in_t, within=c("trial_type", "phase"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values

m1 <- aov(duration ~ (trial_type*phase) + Error(ID/(trial_type*phase)), mean_ex_in_t)
summary(m1)
TukeyHSD(m1)

mean_ex_in_t %>%
  ggplot(aes(y=duration, x=trial_type, fill=phase)) + 
  geom_boxplot() +
  ggtitle("Respiration phase duration") +
  ylab("t in s") +
  #xlab("Time since R-peak in ms")
  theme_cowplot() 

