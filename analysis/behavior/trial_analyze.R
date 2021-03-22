function (d) {
# input "d" <- output by "trial_filter.R"

# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      May 28, 2020


# BLOCK FILTER ------------------------------------------------------------

# Minimum percentage points of hit rate (HR) > false alarm rate (FAR)
margin_hit_FA = 0.05

# Filter for hit rate > false alarm rate
for (i in unique(d$ID)) {
  for (b in unique(d$block[d$ID == i])) {
    HR <- mean(d$resp1[d$ID == i & d$block == b & d$resp_filter == 1 & d$stim_type == 1])
    FAR <- mean(d$resp1[d$ID == i & d$block == b & d$resp_filter == 1 & d$stim_type == 0])
    
    d$HR_larger_FAR_filter[d$ID == i & d$block == b] <- as.integer((HR-FAR)>margin_hit_FA)
  }
}


# CONDTION LABELS --------------------------------------------------------

# Codes to select trials
# Stimulus conditions (stim_type)
null_code <- 0
near_code <- 1

# Response conditions (resp1 & resp2)
no_code <- 0
yes_code <- 1
unconf_code <- 0
conf_code <- 1


# DETECTION & CONFIDENCE MEASURES ------------------------------------------------------

# Prepare data frame that has size of "results"
# Size = "Last line number" - "First line number" + 1
dID <- data.frame(numeric(0),
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

# Loop participants
for (i in unique(d$ID)) {

    # Get valid trial data of one participant:
    # filter for ID, valid responses (see "trial_filter.R") and valid blocks (see above)
    d_valid <- subset(d, ID==i & resp_filter==1 & HR_larger_FAR_filter==1)
    
    # Logic for variable names <stimulus condition>_<variable>_<sub conditions>
    
    results <- c(ID = i,
                 valid_blocks = as.integer(paste(sort(unique(d_valid$block)), collapse = '')),
                 block_num = length(unique(d_valid$block)),
                 trial_t = mean(d_valid$trial_t, na.rm = TRUE),
                 ISI = mean(d_valid$ISI, na.rm = TRUE),
                 null_trial_t = mean(d_valid$trial_t[d_valid$stim_type == null_code], na.rm = TRUE),
                 null_trial_t_no = mean(d_valid$trial_t[d_valid$stim_type == null_code & d_valid$resp1 == no_code], na.rm = TRUE),
                 null_trial_t_yes = mean(d_valid$trial_t[d_valid$stim_type == null_code & d_valid$resp1 == yes_code], na.rm = TRUE),
                 null_intensity = mean(d_valid$intensity[d_valid$stim_type == null_code]),
                 null_resp1 = mean(d_valid$resp1[d_valid$stim_type == null_code]),
                 null_resp1_unconf = mean(d_valid$resp1[d_valid$stim_type == null_code & d_valid$resp2 == unconf_code]),
                 null_resp1_conf = mean(d_valid$resp1[d_valid$stim_type == null_code & d_valid$resp2 == conf_code]),
                 null_resp2 = mean(d_valid$resp2[d_valid$stim_type == null_code]),
                 null_resp2_no = mean(d_valid$resp2[d_valid$stim_type == null_code & d_valid$resp1 == no_code]),
                 null_resp2_yes = mean(d_valid$resp2[d_valid$stim_type == null_code & d_valid$resp1 == yes_code]),
                 
                 null_resp1_t = mean(d_valid$resp1_t[d_valid$stim_type == null_code]),
                 null_resp1_t_yes = mean(d_valid$resp1_t[d_valid$stim_type == null_code & d_valid$resp1 == yes_code]),
                 null_resp1_t_no = mean(d_valid$resp1_t[d_valid$stim_type == null_code & d_valid$resp1 == no_code]),
                 null_resp1_t_conf = mean(d_valid$resp1_t[d_valid$stim_type == null_code & d_valid$resp2 == conf_code]),
                 null_resp1_t_unconf = mean(d_valid$resp1_t[d_valid$stim_type == null_code & d_valid$resp2 == unconf_code]),
                 null_resp1_t_no_conf = mean(d_valid$resp1_t[d_valid$stim_type == null_code& d_valid$resp1 == no_code & d_valid$resp2 == conf_code]),
                 null_resp1_t_no_unconf = mean(d_valid$resp1_t[d_valid$stim_type == null_code& d_valid$resp1 == no_code & d_valid$resp2 == unconf_code]),
                 
                 null_resp2_t = mean(d_valid$resp2_t[d_valid$stim_type == null_code]),
                 null_resp2_t_no = mean(d_valid$resp2_t[d_valid$stim_type == null_code & d_valid$resp1 == no_code]),
                 null_resp2_t_conf = mean(d_valid$resp2_t[d_valid$stim_type == null_code & d_valid$resp2 == conf_code]),
                 null_resp2_t_unconf = mean(d_valid$resp2_t[d_valid$stim_type == null_code & d_valid$resp2 == unconf_code]),
                 null_resp2_t_no_conf = mean(d_valid$resp2_t[d_valid$stim_type == null_code& d_valid$resp1 == no_code & d_valid$resp2 == conf_code]),
                 null_resp2_t_no_unconf = mean(d_valid$resp2_t[d_valid$stim_type == null_code& d_valid$resp1 == no_code & d_valid$resp2 == unconf_code]),
                 
                 null_n = sum(d_valid$stim_type == null_code), #, na.rm = TRUE
                 null_n_no = sum(d_valid$stim_type == null_code & d_valid$resp1 == no_code),
                 null_n_yes = sum(d_valid$stim_type == null_code & d_valid$resp1 == yes_code),
                 null_n_conf = sum(d_valid$stim_type == null_code & d_valid$resp2 == conf_code),
                 null_n_unconf = sum(d_valid$stim_type == null_code & d_valid$resp2 == unconf_code),
                 null_n_no_conf = sum(d_valid$stim_type == null_code & d_valid$resp1 == no_code & d_valid$resp2 == conf_code),
                 null_n_no_unconf = sum(d_valid$stim_type == null_code & d_valid$resp1 == no_code & d_valid$resp2 == unconf_code),
                 null_n_yes_conf = sum(d_valid$stim_type == null_code & d_valid$resp1 == yes_code & d_valid$resp2 == conf_code),
                 null_n_yes_unconf = sum(d_valid$stim_type == null_code & d_valid$resp1 == yes_code & d_valid$resp2 == unconf_code),
                 
                 near_trial_t = mean(d_valid$trial_t[d_valid$stim_type == near_code], na.rm = TRUE),
                 near_trial_t_no = mean(d_valid$trial_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code], na.rm = TRUE),
                 near_trial_t_yes = mean(d_valid$trial_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code], na.rm = TRUE),
                 near_intensity = mean(d_valid$intensity[d_valid$stim_type == near_code]),
                 near_resp1 = mean(d_valid$resp1[d_valid$stim_type == near_code]),
                 near_resp1_unconf = mean(d_valid$resp1[d_valid$stim_type == near_code & d_valid$resp2 == unconf_code]),
                 near_resp1_conf = mean(d_valid$resp1[d_valid$stim_type == near_code & d_valid$resp2 == conf_code]),
                 near_resp2 = mean(d_valid$resp2[d_valid$stim_type == near_code]),
                 near_resp2_no = mean(d_valid$resp2[d_valid$stim_type == near_code & d_valid$resp1 == no_code]),
                 near_resp2_yes = mean(d_valid$resp2[d_valid$stim_type == near_code & d_valid$resp1 == yes_code]),
                 near_n = sum(d_valid$stim_type == near_code),
                 near_n_yes = sum(d_valid$stim_type == near_code & d_valid$resp1 == yes_code),
                 near_n_no = sum(d_valid$stim_type == near_code & d_valid$resp1 == no_code),
                 near_n_conf = sum(d_valid$stim_type == near_code & d_valid$resp2 == conf_code),
                 near_n_unconf = sum(d_valid$stim_type == near_code & d_valid$resp2 == unconf_code),
                 near_n_yes_conf = sum(d_valid$stim_type == near_code & d_valid$resp1 == yes_code & d_valid$resp2 == conf_code),
                 near_n_yes_unconf = sum(d_valid$stim_type == near_code & d_valid$resp1 == yes_code & d_valid$resp2 == unconf_code),
                 near_n_no_conf = sum(d_valid$stim_type == near_code & d_valid$resp1 == no_code & d_valid$resp2 == conf_code),
                 near_n_no_unconf = sum(d_valid$stim_type == near_code & d_valid$resp1 == no_code & d_valid$resp2 == unconf_code),
                 
                 near_trial_median_yes = median(d_valid$trial[d_valid$stim_type == near_code & d_valid$resp1 == yes_code]),
                 
                 near_resp1_t = mean(d_valid$resp1_t[d_valid$stim_type == near_code]),
                 near_resp1_t_no = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code]),
                 near_resp1_t_no_unconf = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code & d_valid$resp2 == unconf_code]),
                 near_resp1_t_no_conf = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code & d_valid$resp2 == conf_code]),
                 near_resp1_t_yes = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code]),
                 near_resp1_t_yes_unconf = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code & d_valid$resp2 == unconf_code]),
                 near_resp1_t_yes_conf = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code & d_valid$resp2 == conf_code]),
                 near_resp1_t_conf = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp2 == conf_code]),
                 near_resp1_t_unconf = mean(d_valid$resp1_t[d_valid$stim_type == near_code & d_valid$resp2 == unconf_code]),
                 
                 near_resp2_t = mean(d_valid$resp2_t[d_valid$stim_type == near_code]),
                 near_resp2_t_no = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code]),
                 near_resp2_t_no_unconf = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code & d_valid$resp2 == unconf_code]),
                 near_resp2_t_no_conf = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp1 == no_code & d_valid$resp2 == conf_code]),
                 near_resp2_t_yes = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code]),
                 near_resp2_t_yes_unconf = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code & d_valid$resp2 == unconf_code]),
                 near_resp2_t_yes_conf = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp1 == yes_code & d_valid$resp2 == conf_code]),
                 near_resp2_t_conf = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp2 == conf_code]),
                 near_resp2_t_unconf = mean(d_valid$resp2_t[d_valid$stim_type == near_code & d_valid$resp2 == unconf_code]))
    
    
    # Update data frame
    dID[nrow(dID)+1,] <- results
    
    # If first row, update column names of data frame
    if (nrow(dID)==1) { colnames(dID) <- names(results)}
    
}

return(dID)

}
