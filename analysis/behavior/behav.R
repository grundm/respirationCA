# respirationCA - Behavioral Data -----------------------------------------------

# Author:         Martin Grund
# Last update:    June 3, 2020

# Mac OS: "alt + cmd + T" runs section
# Linux: "ctrl + alt + T" runs section

# Prepare data (concatenate all txt files of exp_save)
#cat ID*/respirationCA_*_trials_*.txt > respirationCA_trials_tmp.txt
#sed '2,${/ID*/d}' respirationCA_trials_tmp.txt > respirationCA_trials.txt

rm(list = ls())

# SETTINGS ----------------------------------------------------------------

code_dir <- 'respirationca/behavior'
data_dir <- 'respirationca/behavior'

behav_data <- paste(data_dir, '/respirationCA_trials.txt', sep = '')

# PRE-PROCESSING ----------------------------------------------------------

# Load functions
setwd(code_dir)

trial_filter <- dget("trial_filter.R")
block_analyze <- dget("block_analyze.R")
trial_analyze <- dget("trial_analyze.R")
ID_level <- dget("ID_level.R")

# Load data
d <- read.table(behav_data, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

# Trial level
dt <- trial_filter(d)
  
# Block level
db <- block_analyze(dt)

# Participant level
dID1 <- ID_level(db)

# Participant level across blocks
dID <- trial_analyze(dt)

# Demography --------------------------------------------------------------

# Gender
sum(dt$gender[dt$block == 1 & dt$trial == 1]=='w')

# Age
mean(dt$age[dt$block == 1 & dt$trial == 1])
range(dt$age[dt$block == 1 & dt$trial == 1])


# Number of trials --------------------------------------------------------

# Catch trials
mean(dID$null_n)
range(dID$null_n)

# Near-threshold trials
mean(dID$near_n)
range(dID$near_n)

# Correct rejections
mean(dID$null_n_no)
range(dID$null_n_no)

# Hits
mean(dID$near_n_yes)
range(dID$near_n_yes)

# Misses
mean(dID$near_n_no)
range(dID$near_n_no)

# False alarms
mean(dID$null_n_yes)
range(dID$null_n_yes)
sum(dID$null_n_yes == 0) # Number of participants without false alarms

# Confident misses
mean(dID$near_n_no_conf)
range(dID$near_n_no_conf)

# Unconfident misses
mean(dID$near_n_no_unconf)
range(dID$near_n_no_unconf)

# Confident hits
mean(dID$near_n_yes_conf)
sd(dID$near_n_yes_conf)
range(dID$near_n_yes_conf)

# Unconfident hits
mean(dID$near_n_yes_unconf)
sd(dID$near_n_yes_unconf)
range(dID$near_n_yes_unconf)

t.test(dID$near_n_yes_conf,dID$near_n_yes_unconf,paired = T)

# Confident correct rejections
mean(dID$null_n_no_conf)
range(dID$null_n_no_conf)

# Unconfident correct rejections
mean(dID$null_n_no_unconf)
range(dID$null_n_no_unconf)
sum(dID$null_n_no_unconf == 0) # 2

# Unconfident false alarms
mean(dID$null_n_yes_unconf)
range(dID$null_n_yes_unconf)
sum(dID$null_n_yes_unconf == 0) # 6

# Confident false alarms
mean(dID$null_n_yes_conf)
sd(dID$null_n_yes_conf)
range(dID$null_n_yes_conf)
sum(dID$null_n_yes_conf == 0) # 12


# Intensity ---------------------------------------------------------------

mean(dID$near_intensity)
range(dID$near_intensity)

# participants with SEP recording
mean(dID$near_intensity[dID$ID %in% c(6:17)])
range(dID$near_intensity[dID$ID %in% c(6:17)])

# Trial length ------------------------------------------------------------

# All
mean(dID$trial_t)
range(dID$trial_t)

# Difference between trial length and ISI?
mean(dID$ISI)
