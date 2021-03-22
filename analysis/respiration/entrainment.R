###########################################################################
####                            Entrainment                            ####
###########################################################################

# Author:         Martin Grund, Marc Pabst
# Last update:    February 10, 2021


# Settings ---------------------------------------------------------------

rm(list=ls())

data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

library(circular)

# Load data ---------------------------------------------------------------

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210208.Rds", sep="/"))


# Select valid data -------------------------------------------------------

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_filter == 1)


# Difference to block mean ------------------------------------------------

data_stims_valid$stim_degree <- circular(data_stims_valid$stim_degree, type="angles", units="degree", modulo="2pi", rotation="clock", zero=0)

# Calculate mean angle in first block for each participant
block_angle_mean <- aggregate(stim_degree ~ ID, subset(data_stims_valid, block == 1), FUN = mean)

# Define as circular object
block_angle_mean$stim_degree <- circular(block_angle_mean$stim_degree, type="angles", units="degree", modulo = "2pi", rotation="clock", zero=0)

# Prepare difference to mean calculation
data_stims_valid$diff_angle2mean <- NA

for (i in block_angle_mean$ID) {
  print(i)
  
  diff <- data_stims_valid$stim_degree[data_stims_valid$ID == i & data_stims_valid$block == 1] - block_angle_mean$stim_degree[block_angle_mean$ID == i]
  
  diff2 <- abs(diff) %% 360

  diff3 <- ifelse(diff2 > 180, 360 - diff2, diff2)
  
  data_stims_valid$diff_angle2mean[data_stims_valid$ID == i & data_stims_valid$block == 1] <- diff3
  
}


# Model -------------------------------------------------------------------

# Load packages
library("tidyverse")
library("rstatix")
library("cowplot")

# Create Entrainment Data
entrainment_data = data_stims_valid %>% 
  filter(trial_type != "FA") %>%
  filter(block == 1) %>%
  filter(trial < 31) # %>% filter(trial > 1)

trial_num_per_trial <- aggregate(block ~ trial, entrainment_data, FUN="sum")

trial_num_per_ID <- aggregate(block ~ ID, entrainment_data, FUN="sum")

# Excluded outliers with very few trials
# is_outlier requires library("rstatix")
sort(trial_num_per_ID$block)
outlier_ID <- trial_num_per_ID$ID[is_outlier(trial_num_per_ID$block)]

# Exclude only ID22, because 16 trials of ID39 spread across 30 trials
outlier_ID <- 22

`%notin%` <- Negate(`%in%`)

entrainment_data = entrainment_data %>%
  filter(ID %notin% outlier_ID)

# Intercept-only Model
model1 = lme4::lmer(formula = diff_angle2mean ~ 1 + (1 | ID), data = entrainment_data, REML = FALSE)
# Saturated Model
model2 = lme4::lmer(formula = diff_angle2mean ~ 1 + trial + (1 | ID), data = entrainment_data, REML = FALSE)

coef_sum <- coef(summary(model2))
coefs <- coef(model2)
intercepts <- coefs[[1]]$`(Intercept)`

sd(intercepts)/sqrt(length(intercepts))

# Compare Models
anova(model1, model2)

# Plot Entrainment Data 
ggplot(entrainment_data, aes(y = diff_angle2mean, x = trial)) + 
  #  geom_smooth(aes(color = ID), alpha=0.5, se=F, method = "lm",  size = 0.4) +
  geom_boxplot(aes(group = trial, color = ID), fill = "#EB726D", alpha=1) +
  theme_cowplot() + 
  theme(legend.position = "none") +
  geom_abline(slope = coef_sum[2,1], intercept = coef_sum[1,1], lwd = 1.5, alpha = 0.7) +
  theme(axis.title.y = element_blank()) +
  xlab("Trial") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Difference to mean angle per trial")
