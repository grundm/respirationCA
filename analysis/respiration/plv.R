###########################################################################
####                    Phase Locking Value                            ####
###########################################################################

# Author:         Martin Grund
# Last update:    February 11, 2021


# Settings ---------------------------------------------------------------

rm(list=ls())

code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

library(dplyr)
library(rstatix)
library(magrittr)
library(circular)
library("ggpubr")

# Create colormap
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=4
colmap <- c(rgb(141,211,199, maxColorValue=255),
            rgb(255,255,179, maxColorValue=255),
            rgb(190,186,218, maxColorValue=255),
            rgb(251,128,114, maxColorValue=255))


# Load data ---------------------------------------------------------------

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210208.Rds", sep="/"))

ecg_data <- read.csv2(paste(data_dir, 'ecgdata_all.csv', sep = '/'))


# Filter respiratory data -------------------------------------------------

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_filter == 1)


# Match respiratory with ECG data -----------------------------------------------------

resp_data_valid <- data_stims_valid

resp_data_valid$ecg_degree <- NA

resp_data_valid$RR_int_zero_stimulus <- NA

for (i in 1:nrow(resp_data_valid)) {
  
  # Check if corresponding trial in ECG data exists to add degree in cardiac cycle, if not set NA
  resp_data_valid$ecg_degree[i] <- ifelse(length(ecg_data$stim_degree[ecg_data$ID == resp_data_valid$ID[i]
                                                                      & ecg_data$block == resp_data_valid$block[i]
                                                                      & ecg_data$trial == resp_data_valid$trial[i]]),
                                          ecg_data$stim_degree[ecg_data$ID == resp_data_valid$ID[i]
                                                               & ecg_data$block == resp_data_valid$block[i]
                                                               & ecg_data$trial == resp_data_valid$trial[i]],
                                          NA)
  
  # Check if corresponding trial in ECG data exists to add degree in cardiac cycle, if not set NA
  resp_data_valid$RR_int_zero_stimulus[i] <- ifelse(length(ecg_data$RR_int_zero_stimulus[ecg_data$ID == resp_data_valid$ID[i]
                                                                                         & ecg_data$block == resp_data_valid$block[i]
                                                                                         & ecg_data$trial == resp_data_valid$trial[i]]),
                                                    ecg_data$RR_int_zero_stimulus[ecg_data$ID == resp_data_valid$ID[i]
                                                                                  & ecg_data$block == resp_data_valid$block[i]
                                                                                  & ecg_data$trial == resp_data_valid$trial[i]],
                                                    NA)

}


# Estimate respiratory and cardiac frequency ------------------------------

# Use respiratory cycle duration associated with each trial
resp_data_valid$resp_rate <- 1/resp_data_valid$resp_cycle_t

# Mean for each participant across experiment
mean_resp_rate <- aggregate(resp_rate ~ ID, resp_data_valid, FUN=mean)

# Same for cardiac cycle
resp_data_valid$cardiac_rate <- 1/resp_data_valid$RR_int_zero_stimulus

# Mean for each participant across experiment
mean_cardiac_rate <- aggregate(cardiac_rate ~ ID, resp_data_valid, FUN=mean)

# Calculate ratio between mean respiratory and cardiac frequency
m <- mean_resp_rate
m$mean_cardiac <- mean_cardiac_rate$cardiac_rate

m$m <- round(m$mean_cardiac/m$resp_rate)


# Weighted phase difference (respiration - cardiac) -----------------------

resp_data_valid$delta_phase <- NA

for (i in unique(resp_data_valid$ID)) {
  
  resp_data_valid$delta_phase[resp_data_valid$ID == i] <- m$m[m$ID == i] * resp_data_valid$stim_degree[resp_data_valid$ID == i] * (pi/180) - resp_data_valid$ecg_degree[resp_data_valid$ID == i] * (pi/180)

}


# Calculate PLV -----------------------------------------------------------

# Create exponential
resp_data_valid$delta_phase_cpl <- exp(complex(real = 0, imaginary = resp_data_valid$delta_phase))

# Average for each participant and trial type
delta_phase_cpl_mean <- aggregate(delta_phase_cpl ~ ID*trial_type, resp_data_valid, FUN=mean)

# Take absolute -> PLV
delta_phase_cpl_mean$PLV <- abs(delta_phase_cpl_mean$delta_phase_cpl)


# Plot and test PLV -------------------------------------------------------

# Sort 
delta_phase_cpl_mean <- subset(delta_phase_cpl_mean, trial_type != 'FA')

delta_phase_cpl_mean$trial_type <- factor(delta_phase_cpl_mean$trial_type,
                                          levels = c('CR', 'near_miss', 'near_hit'),
                                          ordered = T)

fit <- aov_ez("ID","PLV", delta_phase_cpl_mean, within=c("trial_type"))
fit


# T-test PLVs between conditions
stat_test = delta_phase_cpl_mean %>% 
  t_test(PLV ~ trial_type, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>%
  add_xy_position(step.increase=.05)


ggplot(delta_phase_cpl_mean,aes(x=trial_type, y=PLV)) + 
  geom_boxplot(aes(x=trial_type, y=PLV,  fill=trial_type )) + 
  stat_pvalue_manual(stat_test, label = "p.adj", tip.length = 0.01, hide.ns=T) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()) +
  ylab("Phase-locking value") +
  ggtitle("Cardiac and respiratory coupling") +
  scale_x_discrete(labels = c('CR', 'Miss', 'Hit')) +
  scale_fill_manual(values=colmap[c(1,3,4)])
