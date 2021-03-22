# respirationCA - Pre-/post-stimulus cardiac interval data for hits and misses -----

# Author:         Martin Grund & Esra Al (based on Motyka et al., 2019, https://github.com/Pawel-Motyka/CCSomato/blob/master/_scripts/CCSomato_analysis_main.Rmd)
# Last update:    March 11, 2021



# Settings ----------------------------------------------------------------

rm(list = ls())

# Data path
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

# Code dir
code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'

# Input data (output from ecg_stim2peak.R)
data_file <- paste(data_dir, '/ecgdata_all.csv', sep = '')

# load required packages
require(ggplot2, warn.conflicts = FALSE, quietly=TRUE)
require(scales, warn.conflicts = FALSE, quietly=TRUE)
require(lme4, warn.conflicts = FALSE, quietly=TRUE)
require(tidyr, warn.conflicts = FALSE, quietly=TRUE)
require(afex, warn.conflicts = FALSE, quietly=TRUE)
require(plyr, warn.conflicts = FALSE, quietly=TRUE)

# Load data
data <- read.csv2(data_file)

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

# Minimum percentage points of hit rate (HR) > false alarm rate (FAR)
margin_hit_FA = 0.05


# FILTER ----------------------------------------------------

# Response time filter
data$resp_t_filter <- as.integer(data$resp1_t > resp1_min & data$resp1_t < resp1_max & data$resp2_t > resp2_min & data$resp2_t < resp2_max)

# Response button filter
data$resp_btn_filter <- as.integer(match(data$resp1_btn, resp1_btn, 0) & match(data$resp2_btn, resp2_btn, 0))

# Combine response time & button filter
data$resp_filter <- as.integer(data$resp_t_filter & data$resp_btn_filter)

# Cue presentation time filter
data$cue_t_filter <- as.integer( (data$onset_resp1 - data$onset_cue) > cue_t_min & (data$onset_resp1 - data$onset_cue) < cue_t_max)

# Combine response time & button filter
data$resp_filter <- as.integer(data$resp_t_filter & data$resp_btn_filter)


# Filter for hit rate > false alarm rate
for (i in unique(data$ID)) {
  for (b in unique(data$block[data$ID == i])) {
    HR <- mean(data$resp1[data$ID == i & data$block == b & data$resp_filter == 1 & data$stim_type == 1])
    FAR <- mean(data$resp1[data$ID == i & data$block == b & data$resp_filter == 1 & data$stim_type == 0])
    
    data$HR_larger_FAR_filter[data$ID == i & data$block == b] <- as.integer((HR-FAR)>margin_hit_FA)
  }
}


# Interbeat intervals hit/miss/CR -----------------------------------------

# Specify a data frame with pre- and post-stimulus RR intervals
data_ID <- data.frame(numeric(0),
                      numeric(0), # S-2 hits
                      numeric(0), # S-1 hits
                      numeric(0), # Stimulus hits
                      numeric(0), # S+1 hits
                      numeric(0), # S+2 hits
                      numeric(0), #S-2 misses
                      numeric(0), #S-1 misses
                      numeric(0), #Stimulus misses
                      numeric(0), #S+1 misses
                      numeric(0), #S+2 misses
                      numeric(0), #S-2 CR
                      numeric(0), #S-1 CR
                      numeric(0), #Stimulus CR
                      numeric(0), #S+1 CR
                      numeric(0)) #S+2 CR


for ( i in unique(data$ID )) { 
  
  # Calculate RR intervals
  results <- c(ID = i,
               T0_hit = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_hit = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp1 == 1  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_hit = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp1 == 1  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_hit = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp1 == 1  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_hit = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp1 == 1  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               
               T0_miss = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp1 == 0  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_miss = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp1 == 0  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_miss = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp1 == 0  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_miss = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp1 == 0  & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_miss = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               
               T0_CR = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 0 & data$resp1 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_CR = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 0 & data$resp1 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_CR = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 0 & data$resp1 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_CR = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 0 & data$resp1 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_CR = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 0 & data$resp1 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T))
  
  # Update data frame
  data_ID[nrow(data_ID)+1,] <- results
  
  # If first row, update column names of data frame
  if (nrow(data_ID)==1) { colnames(data_ID) <- names(results)} 
  
}



# Interbeat intervals confident/unconfident misses/hits -------------------
# (Confidence x awareness effect)

# Specify a data frame with pre- and post-stimulus RR intervals
data_ID <- data.frame(numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0))


for ( i in unique(data$ID )) { 
  
  # Calculate RR intervals
  results <- c(ID = i,
               T0_hit_conf = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_hit_conf = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_hit_conf = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_hit_conf = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_hit_conf = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               
               T0_hit_unconf = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_hit_unconf = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_hit_unconf = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_hit_unconf = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_hit_unconf = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp1 == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               
               T0_miss_conf = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_miss_conf = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_miss_conf = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_miss_conf = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_miss_conf = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               
               T0_miss_unconf = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_miss_unconf = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_miss_unconf = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_miss_unconf = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_miss_unconf = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp1 == 0 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T))
  
  # Update data frame
  data_ID[nrow(data_ID)+1,] <- results
  
  # If first row, update column names of data frame
  if (nrow(data_ID)==1) { colnames(data_ID) <- names(results)} 
  
}


# Interbeat intervals confident/unconfident near-threshold trials ---------

# Specify a data frame with pre- and post-stimulus RR intervals
data_ID <- data.frame(numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0))


for ( i in unique(data$ID )) { 
  
  # Calculate RR intervals
  results <- c(ID = i,
               T0_conf = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_conf = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_conf = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_conf = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_conf = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp2 == 1 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               
               T0_unconf = mean(data$RR_int_minus2[data$ID == i & data$stim_type == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T1_unconf = mean(data$RR_int_minus1[data$ID == i & data$stim_type == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T2_unconf = mean(data$RR_int_zero_stimulus[data$ID == i & data$stim_type == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T3_unconf = mean(data$RR_int_plus1[data$ID == i & data$stim_type == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T),
               T4_unconf = mean(data$RR_int_plus2[data$ID == i & data$stim_type == 1 & data$resp2 == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1], na.rm=T))
  
  # Update data frame
  data_ID[nrow(data_ID)+1,] <- results
  
  # If first row, update column names of data frame
  if (nrow(data_ID)==1) { colnames(data_ID) <- names(results)} 
  
}


# Interbeat intervals confident vs unconfident (w/o false alarms) ---------
# Confidence ONLY effect

# Create filter variable for false alarms
data$FA_filter <- ifelse(data$stim_type==0 & data$resp1==1, 1, 0)

data_select <- subset(data, data$FA_filter == 0 & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1)

# Specify a data frame with pre- and post-stimulus RR intervals
data_ID <- data.frame(numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0),
                      numeric(0))


for ( i in unique(data$ID )) { 
  
  # Calculate RR intervals for each participant
  results <- c(ID = i,
               T0_conf = mean(data$RR_int_minus2[data$ID == i & data$resp2 == 1], na.rm=T),
               T1_conf = mean(data$RR_int_minus1[data$ID == i & data$resp2 == 1], na.rm=T),
               T2_conf = mean(data$RR_int_zero_stimulus[data$ID == i  & data$resp2 == 1], na.rm=T),
               T3_conf = mean(data$RR_int_plus1[data$ID == i & data$resp2 == 1], na.rm=T),
               T4_conf = mean(data$RR_int_plus2[data$ID == i  & data$resp2 == 1], na.rm=T),
               
               T0_unconf = mean(data$RR_int_minus2[data$ID == i & data$resp2 == 0], na.rm=T),
               T1_unconf = mean(data$RR_int_minus1[data$ID == i & data$resp2 == 0], na.rm=T),
               T2_unconf = mean(data$RR_int_zero_stimulus[data$ID == i & data$resp2 == 0], na.rm=T),
               T3_unconf = mean(data$RR_int_plus1[data$ID == i & data$resp2 == 0], na.rm=T),
               T4_unconf = mean(data$RR_int_plus2[data$ID == i & data$resp2 == 0], na.rm=T))
  
  # Update data frame
  data_ID[nrow(data_ID)+1,] <- results
  
  # If first row, update column names of data frame
  if (nrow(data_ID)==1) { colnames(data_ID) <- names(results)} 
  
}


# ANOVA cofindence/detection x time ---------------------------------------

## Prepare the data frame for two-way repeated measures ANOVA

library(tidyr)

# Reshape the data frame
dat <- gather(data_ID, key = "time", value = "value", -ID)

# Create separate factor: awareness (hit/miss)
dat <- separate(dat, time, into = c('time','cond'), sep = '_', extra = 'merge')

# Save ID, time & condition as factors
dat$ID <- factor(dat$ID)
dat$time <- factor(dat$time)
dat$cond <- factor(dat$cond)

## Perform two-way repeated measures ANOVA

# Run to see uncorrected degrees of freedom
fit_all <- aov_ez("ID","value", dat, within=c("time", "cond"), return = "nice")
fit_all 

# Perform test with corrected degrees of freedom
fit_all <- aov_ez("ID","value", dat, within=c("time", "cond"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values


  # fit_all$`Sphericity Corrections`$`p[GG]`
  # fit_all$anova_table$`num Df`
  # fit_all$anova_table$`den Df`
  # fit_all$anova_table$F
  # fit_all$anova_table
  # # But Mauchly's Test for Sphericity significant
  # # Report Greenhouse-Geisser correction
  # # Main effect time: F(2.73,109.36) = 35.60, p < 0.0001 (3 * 10^-15)
  # # Main effect stimulus-response condition: F(1.42,56.75) = 2.35, p = 0.12
  # # Main effect interaction: F(4.89,195.55) = 4.92, p = 0.0003



# ANOVA confidence x detection x time -------------------------------------

## Prepare the data frame for two-way repeated measures ANOVA

# Reshape the data frame
dat <- gather(data_ID, key = "time", value = "value", -ID)

# Create separate factor: awareness (hit/miss)
dat <- separate(dat, time, into = c('time','cond','conf'), sep = '_', extra = 'merge')

# Save ID, time & condition as factors
dat$ID <- factor(dat$ID)
dat$time <- factor(dat$time)
dat$cond <- factor(dat$cond)
dat$conf <- factor(dat$conf)

## Perform two-way repeated measures ANOVA

# Run to see uncorrected degrees of freedom
fit_all <- aov_ez("ID","value", dat, within=c("time", "cond", "conf"), return = "nice")
fit_all 

# Perform test with corrected degrees of freedom
fit_all <- aov_ez("ID","value", dat, within=c("time", "cond", "conf"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values

# 
# fit_all$`Sphericity Corrections`$`p[GG]`
# fit_all$anova_table$`num Df`
# fit_all$anova_table$`den Df`
# fit_all$anova_table$F
# fit_all$anova_table


# Define contrasts for custum post-hoc t-tests ----------------------------

# CR vs. miss vs. hit
cond_contrast = list(c('hit','miss'),
                     c('hit','CR'),
                     c('miss','CR'))
# Conf vs. unconf
cond_contrast = list(c('conf','unconf'))

# Conf miss/hit vs. unconf miss/hit
# Create one conditions factor out of detection and confidence factor
dat <- unite(dat, cond, cond:conf)

cond_contrast = list(c('hit_conf','hit_unconf'),
                     c('hit_conf','miss_conf'),
                     c('hit_conf','miss_unconf'),
                     c('hit_unconf','miss_conf'),
                     c('hit_unconf','miss_unconf'),
                     c('miss_conf','miss_unconf'))

# Run post-hoc t-tests ----------------------------------------------------

# FIRST: Compare conditions pairwise at each time point (number of conditions x 5)

# Prepare data frame for p-values
pval <- data.frame(row.names = unique(dat$time))

# Prepare vector for column names
col_names <- c(0)

# Loop time points
for (t in 1:nrow(pval)) {

  # Loop contrasts
  for (c in 1:length(cond_contrast)) {
    
    # Paired t-test
    result <- t.test(dat$value[dat$time == rownames(pval)[t] & dat$cond == cond_contrast[[c]][1]],
                     dat$value[dat$time == rownames(pval)[t] & dat$cond == cond_contrast[[c]][2]],
                     paired = T)
    
    # Save p-value
    pval[t,c] <- result$p.value
  
    # Create string for compared conditions
    col_names[c] <- paste(cond_contrast[[c]][1], cond_contrast[[c]][2], sep = '_')
  }
  
}

# Replace column names with string for compared conditions
colnames(pval) <- col_names

# Transform to long format
pval_long <- gather(pval, key = "cond", value = "value")


# SECOND: Compare subsequent time points for each condition (4 x number of conditions)

pval2 <- data.frame(row.names = c('T0_T1', 'T1_T2', 'T2_T3', 'T3_T4'))
    
col_names2 <- c(0)

# Loop time points
for (t in 1:nrow(pval2)) {
  
  for (c in 1:length(unique(dat$cond))) {
    result <- t.test(dat$value[dat$time == unique(dat$time)[t] & dat$cond == unique(dat$cond)[c]],
                     dat$value[dat$time == unique(dat$time)[(t+1)] & dat$cond == unique(dat$cond)[c]],
                     paired = T)
    
    pval2[t,c] <- result$p.value
    
    col_names2[c] <- as.character(unique(dat$cond)[c])
  }
  
}

colnames(pval2) <- col_names2

pval2_long <- gather(pval2, key = "cond", value = "value")

# FDR-correction
p_adj <- p.adjust(c(pval_long$value, pval2_long$value), method = "fdr")

# Which tests are significant
which(p_adj < 0.05)


# CR, miss, & hit ---------------------------------------------------------

p_adj[16:27] < 0.005

t.test(dat$value[dat$time == "T3" & dat$cond == "conf"]*1000,
       dat$value[dat$time == "T3" & dat$cond == "unconf"]*1000,
       paired = T)

# T3: hit vs. miss
# T3: hit vs. CR
# T4: hit vs. CR

# hit: all timepoints
# miss: all timepoints
# CR: all timepoints


# Confidence intervals for two-way ANOVA ----------------------------------

## Create three functions used to calculate within-subjects variation. These fragments of code are credited to Winston Chang and were used and copied here under the CC0 license (http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#understanding-within-subjects-error-bars)
# Get within variance estimate helper functions
source(paste(code_dir, 'assets/CI_within_helper.R', sep = '/'))

## Calculate within-subjects variation
d <- summarySEwithin(data=dat, measurevar = 'value', betweenvars=NULL, withinvars=c('time','cond'), idvar='ID', na.rm=FALSE, conf.interval=.95, .drop=TRUE)

# Change into milliseconds 
d$value <- d$value * 1000
d$ci <- d$ci * 1000


# Confidence intervals for three-way ANOVA --------------------------------

## Create three functions used to calculate within-subjects variation. These fragments of code are credited to Winston Chang and were used and copied here under the CC0 license (http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#understanding-within-subjects-error-bars)
# Get within variance estimate helper functions
source(paste(code_dir, 'assets/CI_within_helper.R', sep = '/'))

## Calculate within-subjects variation
d <- summarySEwithin(data=dat, measurevar = 'value', betweenvars=NULL, withinvars=c('time','cond','conf'), idvar='ID', na.rm=FALSE, conf.interval=.95, .drop=TRUE)

# Change into milliseconds 
d$value <- d$value * 1000
d$ci <- d$ci * 1000

# PLOT RR INTERVALS -------------------------------------------------------

# Create colormap
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=4
colmap <- c(rgb(141,211,199, maxColorValue=255),
            rgb(255,255,179, maxColorValue=255),
            rgb(190,186,218, maxColorValue=255),
            rgb(251,128,114, maxColorValue=255))

# https://colorbrewer2.org/?type=qualitative&scheme=Dark2&n=3
colmap <- c(rgb(27,158,119, maxColorValue=255),
            rgb(117,112,179, maxColorValue=255),
            rgb(217,95,2, maxColorValue=255))


#colmap <- c("yellow", "firebrick3", "blue3")

# Prepare descriptions
labs <- c("S-2", "S-1", "Stimulus", "S+1", "S+2")

# Save awareness as factor
d$cond <- factor(d$cond, levels = c('CR','miss','hit'))

## Plot associations between heart rate and perceptual performance over the course of a trial
pd <- position_dodge(0)
plot <- ggplot(d, aes(x=time, y=value,group=cond,colour=cond)) + 
  geom_ribbon(aes(ymin=value-ci, ymax=value+ci, fill = factor(cond)), position=pd, show.legend = F, alpha = 0.4, colour = NA) +
  scale_fill_manual(values = colmap[c(1,2,3)]) +
  geom_line(position=pd, size = 0.7) +
  geom_point(position=pd) + 
  theme_classic() +
  labs(x = " ", y = "Interbeat interval (ms)") +
  scale_x_discrete(labels=labs) +
  scale_y_continuous(limits = c(820,860)) +
  scale_color_manual(values = colmap[c(1,2,3)]) +
  coord_fixed(ratio = 0.05) +
  theme(legend.position = c(0.07,0.87), legend.spacing.y = unit(0.2, 'cm'), legend.key.height = unit(0.5, 'cm')) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.key.width = unit(0.6, "cm")) +
  theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.text=element_text(size = 14, colour = "black")) +
  annotate(geom = "text", x = 4, y = 850.5, label = "**", color = "black", size = 7) +
  annotate(geom = "text", x = 5, y = 841.5, label = "*", color = "black", size = 7) +
  annotate("segment", x = 1, xend = 2, y = 825, yend = 825, colour = "black") +
  annotate(geom = "text", x = 1.5, y = 822, label = "*", color = "black", size = 6) +
  annotate("segment", x = 2, xend = 3, y = 824, yend = 824, colour = "black") +
  annotate(geom = "text", x = 2.5, y = 821, label = "*", color = "black", size = 6) +
  annotate("segment", x = 3, xend = 4, y = 825, yend = 825, colour = "black") +
  annotate(geom = "text", x = 3.5, y = 822, label = "*", color = "black", size = 6) +
  annotate("segment", x = 4, xend = 5, y = 824, yend = 824, colour = "black") +
  annotate(geom = "text", x = 4.5, y = 821, label = "*", color = "black", size = 6)
print(plot)

# PLOT RR INTERVALS - CONFIDENCE X AWARENESS -------------------------------------------------------

# Create colormap
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=4
colmap <- c(rgb(141,211,199, maxColorValue=255),
            rgb(255,255,179, maxColorValue=255),
            rgb(190,186,218, maxColorValue=255),
            rgb(251,128,114, maxColorValue=255))

colmap <- c("firebrick3", "yellow", "blue3", "green3")

# Prepare descriptions
labs <- c("S-2", "S-1", "Stimulus", "S+1", "S+2")

# Save awareness as factor
d$cond <- factor(d$cond)
d$conf <- factor(d$conf)

d$cond2 <- paste(d$cond, d$conf, sep = '_')

## Plot associations between heart rate and perceptual performance over the course of a trial
pd <- position_dodge(0)
plot <- ggplot(d, aes(x=time, y=value,group=cond2, colour=cond2)) + 
  geom_ribbon(aes(ymin=value-ci, ymax=value+ci, fill = factor(cond2)), position=pd, show.legend = F, alpha = 0.3, colour = NA) +
  scale_fill_manual(values = colmap[c(1,2,3,4)]) +
  geom_line(position=pd, size = 0.7) +
  geom_point(position=pd) + 
  theme_classic() +
  labs(x = " ", y = "Interbeat interval (ms)") +
  scale_x_discrete(labels=labs) +
  scale_y_continuous(limits = c(820,860)) +
  scale_color_manual(values = colmap[c(1,2,3,4)]) +
  coord_fixed(ratio = 0.05) +
  theme(legend.position = c(0.11,0.87), legend.spacing.y = unit(0.2, 'cm'), legend.key.height = unit(0.5, 'cm')) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.key.width = unit(0.6, "cm")) +
  theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.text=element_text(size = 14, colour = "black")) #+
  #annotate(geom = "text", x = 5, y = 851, label = "*", color = "black", size = 7) +
  #annotate(geom = "text", x = 5, y = 843.5, label = "*", color = "black", size = 7)
print(plot)

t.test(dat$value[dat$time == "T3" & dat$cond == "hit" & dat$conf == "conf"]*1000,
       dat$value[dat$time == "T3" & dat$cond == "miss" & dat$conf == "conf"]*1000,
       paired = T)
# T3 conf: hit vs miss uncorrected: t = 2.1036, df = 40, p-value = 0.04175
# T4 hit: conf vs. unconf uncorrected: t = -2.37, df = 40, p-value = 0.0227
# T4 miss: conf vs. unconf uncorrected: t = -3.2307, df = 40, p-value = 0.002472


# PLOT RR INTERVALS - CONFIDENCE ONLY EFFECT -------------------------------------------------------

# Create colormap
colmap <- c(rgb(31,120,180, maxColorValue=255), "firebrick3")

# Prepare descriptions
labs <- c("S-2", "S-1", "Stimulus", "S+1", "S+2")

# Save awareness as factor
d$cond <- factor(d$cond)
d$conf <- factor(d$conf)

#d$cond2 <- paste(d$cond, d$conf, sep = '_')

d$cond2 <- d$cond

## Plot associations between heart rate and perceptual performance over the course of a trial
pd <- position_dodge(0)
plot <- ggplot(d, aes(x=time, y=value,group=cond2, colour=cond2)) + 
  geom_ribbon(aes(ymin=value-ci, ymax=value+ci, fill = factor(cond2)), position=pd, show.legend = F, alpha = 0.4, colour = NA) +
  scale_fill_manual(values = colmap[c(1,2,3,4)]) +
  geom_line(position=pd, size = 0.7) +
  geom_point(position=pd) + 
  theme_classic() +
  labs(x = " ", y = "Interbeat interval (ms)") +
  scale_x_discrete(labels=labs) +
  scale_y_continuous(limits = c(820,860)) +
  scale_color_manual(values = colmap[c(1,2,3,4)]) +
  coord_fixed(ratio = 0.05) +
  theme(legend.position = c(0.07,0.87), legend.spacing.y = unit(0.2, 'cm'), legend.key.height = unit(0.5, 'cm')) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.key.width = unit(0.6, "cm")) +
  theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.text=element_text(size = 14, colour = "black")) +
  annotate(geom = "text", x = 4, y = mean(d$value[d$time=="T3"])-1, label = "*", color = "black", size = 7) +
  annotate(geom = "text", x = 5, y = mean(d$value[d$time=="T4"])-1, label = "*", color = "black", size = 7)
print(plot)

t.test(dat$value[dat$time == "T3" & dat$cond == "unconf"]*1000,
       dat$value[dat$time == "T3" & dat$cond == "conf"]*1000,
       paired = T)

t.test(dat$value[dat$time == "T4" & dat$cond == "unconf"]*1000,
       dat$value[dat$time == "T4" & dat$cond == "conf"]*1000,
       paired = T)

