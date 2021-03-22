# respirationCA - ECG plot mean angle for hits and misses -----

# Author:         Martin Grund
# Last update:    August 27, 2020


# Settings ----------------------------------------------------------------

rm(list=ls())

# Input
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'
code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'

# input_file <- paste(data_dir, '/ecg_final.csv', sep = '')
# -> Don't exclude participants with non-uniform distribution of stimulus onsets
input_file <- paste(data_dir, '/ecgdata_all.csv', sep = '')

# Output
fig_path <- '~/ownCloud/promotion/experiment/respirationCA/fig'

circular_pdf <- paste(fig_path, '/circular.pdf', sep = '')

require(circular)
library(ggplot2)


# Load data
data <- read.csv2(input_file)


# FILTER SETTINGS ---------------------------------------------------------

# Response times
resp1_min <- 0 #0.1
resp1_max <- 1.5
resp2_min <- 0 #0.02 # Participants are super fast for resp2
resp2_max <- resp1_max

# Response buttons
# if not one of the buttons pressed, then maximum response time
# sum((d$resp1_btn<1) - (d$resp1_t>=1.5)) -> 0
resp1_btn <- c(1,2)
resp2_btn <- c(1,2)

# Cue presentation time interval
cue_t_min <- 0.61
cue_t_max <- 0.62

# Minimum percentage points of hit rate (HR) > false alarm rate (FAR)
margin_hit_FA = 0.05

# FILTER ----------------------------------------------------

# Response time filter
data$resp_t_filter <- as.integer(data$resp1_t > resp1_min & data$resp1_t < resp1_max & data$resp2_t > resp2_min & data$resp2_t < resp2_max)

# Response button filter
data$resp_btn_filter <- as.integer(match(data$resp1_btn, resp1_btn, 0) & match(data$resp2_btn, resp2_btn, 0))

# Cue presentation time filter
data$cue_t_filter <- as.integer( (data$onset_resp1 - data$onset_cue) > cue_t_min & (data$onset_resp1 - data$onset_cue) < cue_t_max)

# Combine response time & button filter
#data$resp_filter <- as.integer(data$resp_t_filter & data$resp_btn_filter & data$cue_t_filter)
data$resp_filter <- as.integer(data$resp_t_filter & data$resp_btn_filter)

# Filter for hit rate > false alarm rate
# Since trials without ECG recording have been removed in a previous step, there is the possibility of removing blocks
# that would be not removed based on the behavioral data
for (i in unique(data$ID)) {
  for (b in unique(data$block[data$ID == i])) {
    HR <- mean(data$resp1[data$ID == i & data$block == b & data$resp_filter == 1 & data$stim_type == 1])
    FAR <- mean(data$resp1[data$ID == i & data$block == b & data$resp_filter == 1 & data$stim_type == 0])
    
    data$HR_larger_FAR_filter[data$ID == i & data$block == b] <- as.integer((HR-FAR)>margin_hit_FA)
  }
}

# Create data frame with filter
data_valid <- subset(data, resp_filter == 1 & HR_larger_FAR_filter == 1)

# Check all participants --------------------------------------------------

source(paste(code_dir, '/ecg/plot_test_R.R', sep = ''))

# Plot each stimulus-response condition -----------------------------------

pdf('~/ownCloud/promotion/experiment/respirationCA/fig/circular.pdf', width=6, height=6, useDingbats=FALSE)

# Near hit (R = 0.34, p = 0.007)
plot_test_R(data_valid,1,1,c(0,1),'Hit')

# Near miss (R = 0.20, p = 0.18)
plot_test_R(data_valid,1,0,c(0,1),'Miss')

# CR (R = 0.05, p = 0.91)
plot_test_R(data,0,0,c(0,1),'Correct rejection')

# FA (R = 0.34, p = 0.01)
# Excludes 4 participants
plot_test_R(data_valid,0,1,c(0,1),'False alarm')

# Near confident (R = 0.30, p = 0.04)
plot_test_R(data_valid,1,c(0,1),1,'Near confident')

# Near unconfident (R = 0.21, p = 0.20)
plot_test_R(data_valid,1,c(0,1),0,'Near unconfident')

# Near confident miss (R = 0.28, p = 0.04)
plot_test_R(data_valid,1,0,1,'Confident miss')

# Near confident hit (R = 0.38, p = 0.002)
plot_test_R(data_valid,1,1,1,'Confident hit')

# Near unconfident miss (R = 0.22, p = 0.14)
plot_test_R(data_valid,1,0,0,'Unconfident miss')

# Near unconfident hit (R = 0.10, p = 0.69)
plot_test_R(data_valid,1,1,0,'Unconfident hit')

# Confident correct rejection (R = 0.005, p = 1.00)
plot_test_R(data_valid,0,0,1,'Confident correct rejection')

# Unconfident correct rejection (R = 0.01, p = 1.00)
plot_test_R(data_valid,0,0,0,'Unonfident correct rejection')

# Confident false alarms (R = 0.32, p = 0.047)
# Excludes 12 participants
plot_test_R(data_valid,0,1,1,'Confident false alarm')

# Unconfident false alarms (R = 0.17, p = 0.36)
# Excludes 6 participants
plot_test_R(data_valid,0,1,0,'Unconfident false alarms')

dev.off()
