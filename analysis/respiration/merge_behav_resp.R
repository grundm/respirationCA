###########################################################################
####                 Merge Behavior & Respiratory Data                 ####
###########################################################################

# Author:         Martin Grund
# Last update:    June 22, 2021


# Settings ---------------------------------------------------------------

rm(list=ls())

data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'
code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'


# Get data ----------------------------------------------------------------

# Load behavioral data
d_tmp <- read.table(paste(data_dir, 'respirationCA_trials.txt', sep = '/'),
                header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

# Preprocess behavioral data
trial_filter <- dget(paste(code_dir, 'behavior/trial_filter.R', sep = '/'))

data_stims <- trial_filter(d_tmp)

# Load respiratory data
resp_stim_circ = read.csv(paste(data_dir, "resp_stim_circ_20210622.csv", sep="/"),
                          header = F,
                          col.names = c('ID','block','trial',
                                        'diff2inhale_onset','resp_cycle_t_inhale','stim_degree_inhale','resp_cycle_t_next_inhale',
                                        'diff2exhale_onset','resp_cycle_t_exhale','stim_degree_exhale','resp_cycle_t_next_exhale',
                                        'inhale_duration_inhale', 'exhale_duration_inhale',
                                        'inhale_duration_exhale', 'exhale_duration_exhale',
                                        'resp_cycle_t_median_inhale','resp_cycle_t_median_exhale'))


# Merge respiration with behavioral data ----------------------------------

data_stims$diff2inhale_onset <- NA
data_stims$resp_cycle_t_inhale <- NA
data_stims$stim_degree_inhale <- NA
data_stims$resp_cycle_t_next_inhale <- NA

data_stims$diff2exhale_onset <- NA
data_stims$resp_cycle_t_exhale <- NA
data_stims$stim_degree_exhale <- NA
data_stims$resp_cycle_t_next_exhale <- NA

data_stims$inhale_duration_inhale <- NA
data_stims$exhale_duration_inhale <- NA

data_stims$inhale_duration_exhale <- NA
data_stims$exhale_duration_exhale <- NA

data_stims$resp_cycle_t_median_inhale <- NA
data_stims$resp_cycle_t_median_exhale <- NA

for (i in 1:nrow(data_stims)) {
  
  print(i)
  data_stims[i,(ncol(data_stims)-13):ncol(data_stims)] <- resp_stim_circ[resp_stim_circ$ID == data_stims$ID[i]
                                                                        & resp_stim_circ$block == data_stims$block[i]
                                                                        & resp_stim_circ$trial == data_stims$trial[i],
                                                                        4:ncol(resp_stim_circ)]
}


# Filter resp_cycle_median ------------------------------------------------

factor_median_outlier <- 2

#data_stims$resp_cycle_filter <- as.integer(data_stims$resp_cycle_t < data_stims$resp_cycle_t_median*factor_median_outlier)

data_stims$resp_cycle_inhale_filter <- as.integer(data_stims$resp_cycle_t_inhale < data_stims$resp_cycle_t_median_inhale*factor_median_outlier)
data_stims$resp_cycle_exhale_filter <- as.integer(data_stims$resp_cycle_t_exhale < data_stims$resp_cycle_t_median_exhale*factor_median_outlier)


# Add block filter HR > FAR -----------------------------------------------

# Minimum percentage points of hit rate (HR) > false alarm rate (FAR)
margin_hit_FA = 0.05

data_stims$HR_larger_FAR_filter <- NA

# Filter for hit rate > false alarm rate
for (i in unique(data_stims$ID)) {
  for (b in unique(data_stims$block[data_stims$ID == i])) {
    HR <- mean(data_stims$resp1[data_stims$ID == i & data_stims$block == b & data_stims$resp_filter == 1 & data_stims$stim_type == 1])
    FAR <- mean(data_stims$resp1[data_stims$ID == i & data_stims$block == b & data_stims$resp_filter == 1 & data_stims$stim_type == 0])
    
    data_stims$HR_larger_FAR_filter[data_stims$ID == i & data_stims$block == b] <- as.integer((HR-FAR)>margin_hit_FA)
  }
}


# Trial type variable -----------------------------------------------------

data_stims$trial_type <- NA

data_stims$trial_type[data_stims$stim_type==0 & data_stims$resp1==0] <- 'CR'
data_stims$trial_type[data_stims$stim_type==0 & data_stims$resp1==1] <- 'FA'
data_stims$trial_type[data_stims$stim_type==1 & data_stims$resp1==0] <- 'near_miss'
data_stims$trial_type[data_stims$stim_type==1 & data_stims$resp1==1] <- 'near_hit'


# Save data ---------------------------------------------------------------

saveRDS(data_stims,paste(data_dir, "resp_data_stims_martin_20210622.Rds", sep="/"))

