###########################################################################
####                    Phase Locking Value                            ####
###########################################################################

# Author:         Martin Grund
# Last update:    August 6, 2021


# Settings ---------------------------------------------------------------

rm(list=ls())

code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

library(dplyr)
library(rstatix)
library(magrittr)
library(circular)
library(ggpubr)

# Create colormap
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=4
colmap <- c(rgb(141,211,199, maxColorValue=255),
            rgb(255,255,179, maxColorValue=255),
            rgb(190,186,218, maxColorValue=255),
            rgb(251,128,114, maxColorValue=255))


# Load data ---------------------------------------------------------------

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210622.Rds", sep="/")) # Locked to inspiration and expiration

ecg_data <- read.csv2(paste(data_dir, 'ecgdata_all.csv', sep = '/'))


# Filter respiratory data -------------------------------------------------

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_exhale_filter == 1)
data_stims_valid$stim_degree <- data_stims_valid$stim_degree_exhale


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

# IBI across respiration cycle --------------------------------------------

# Distance of stimulus onset to previous expiration onset
resp_data_valid$dist=NaN

# 4 x 1-s interval 0-4 s
intervals <- seq(0,4,1)
# 4 x 90° interval 0-360°
intervals <- seq(0,360,90)
intervals <- seq(0,360,45)

for (i in 1:(length(intervals)-1)) {
  print(i)
  print(intervals[i])
  print(intervals[i+1])
  #data_stims$dist[data_stims$diff2inhale_onset>=intervals[i] & data_stims$diff2inhale_onset<intervals[i+1]] <- i
  #data_stims$dist[data_stims$diff2exhale_onset>=intervals[i] & data_stims$diff2exhale_onset<intervals[i+1]] <- i
  resp_data_valid$dist[resp_data_valid$stim_degree>=intervals[i] & resp_data_valid$stim_degree<intervals[i+1]] <- i
}


# Plot IBI across respiration cycle ---------------------------------------

# Create labels based on intervals
label_list <- paste(as.character(head(intervals,-1)*1000),as.character(intervals[-1]*1000),sep = "\n-")
label_list <- paste(as.character(head(intervals,-1)),'\n-',as.character(intervals[-1]),'°',sep = '')
names(label_list) <- 1:(length(intervals)-1)

mean_IBI <- aggregate(RR_int_zero_stimulus ~ ID*dist, resp_data_valid, FUN='sd')

mean_IBI$dist <- as.factor(mean_IBI$dist)

mean_IBI_base <- mean_IBI
for (i in unique(mean_IBI_base$ID)) {
  mean_IBI_base$RR_int_zero_stimulus[mean_IBI_base$ID == i] <- mean_IBI_base$RR_int_zero_stimulus[mean_IBI_base$ID == i] - mean_IBI_base$RR_int_zero_stimulus[mean_IBI_base$ID == i & mean_IBI_base$dist == 1]
}

stat_test_mean_IBI <- mean_IBI_base %>%
  t_test(RR_int_zero_stimulus ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

# Plot number of false alarms
mean_IBI_base %>%
  ggplot(aes(y=RR_int_zero_stimulus*1000, x=dist)) + 
  geom_boxplot(fill = '#e78ac3') +
  scale_x_discrete(labels=label_list) +
  ggtitle("IBI across respiration cycle") +
  ggtitle("IBI variance across respiration cycle") +
  xlab("Respiration phase") + 
  #ylab('delta IBI in ms to 0-45°') +
  ylab('delta SD in ms to 0-45°') +
  theme_cowplot()
  #stat_pvalue_manual(stat_test_mean_IBI, label = "p.adj", tip.length = 0.01, hide.ns=T)


library("afex")

fit_all <- aov_ez("ID","RR_int_zero_stimulus", mean_IBI, within=c("dist"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values
