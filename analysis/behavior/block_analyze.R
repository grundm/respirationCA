function (d) {
# input "d" <- output by "trial_filter.R"
  
# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      April 28, 2020
  
# FILTER SETTINGS ----------------------------------------------------------------

# Number of trials
null_n_min <-	0 #40
near_n_min <-	0 #90

# Detection
null_resp1_max <- 1#0.20
near_resp1_min <-	0#0.20#0.25
near_resp1_max <- 1#0.80#0.75

margin_hit_FA = 0.05

# Median trial of near hits
near_trial_median_yes_min <- 25#30 # exception ">" not ">=" as usual
near_trial_median_yes_max <- 125#120


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
db <- data.frame(numeric(0),
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
  # Loop blocks
  for (j in unique(d$block) ) {
    
    # Get trial data of one block:
    # filter for ID, block and valid responses (see "trial_filter.R")
    d_block <- subset(d, ID==i & block==j & resp_filter)

    # Logic for variable names <stimulus condition>_<variable>_<sub conditions>
    
    results <- c(ID = i,
                 block = j,
                 trial_t = mean(d_block$trial_t, na.rm = TRUE),
                 ISI = mean(d_block$ISI, na.rm = TRUE),
                 null_trial_t = mean(d_block$trial_t[d_block$stim_type == null_code], na.rm = TRUE),
                 null_trial_t_no = mean(d_block$trial_t[d_block$stim_type == null_code & d_block$resp1 == no_code], na.rm = TRUE),
                 null_trial_t_yes = mean(d_block$trial_t[d_block$stim_type == null_code & d_block$resp1 == yes_code], na.rm = TRUE),
                 null_intensity = mean(d_block$intensity[d_block$stim_type == null_code]),
                 null_resp1 = mean(d_block$resp1[d_block$stim_type == null_code]),
                 null_resp1_unconf = mean(d_block$resp1[d_block$stim_type == null_code & d_block$resp2 == unconf_code]),
                 null_resp1_conf = mean(d_block$resp1[d_block$stim_type == null_code & d_block$resp2 == conf_code]),
                 null_resp2 = mean(d_block$resp2[d_block$stim_type == null_code]),
                 null_resp2_no = mean(d_block$resp2[d_block$stim_type == null_code & d_block$resp1 == no_code]),
                 null_resp2_yes = mean(d_block$resp2[d_block$stim_type == null_code & d_block$resp1 == yes_code]),
                 
                 null_resp1_t = mean(d_block$resp1_t[d_block$stim_type == null_code]),
                 null_resp1_t_yes = mean(d_block$resp1_t[d_block$stim_type == null_code & d_block$resp1 == yes_code]),
                 null_resp1_t_no = mean(d_block$resp1_t[d_block$stim_type == null_code & d_block$resp1 == no_code]),
                 null_resp1_t_conf = mean(d_block$resp1_t[d_block$stim_type == null_code & d_block$resp2 == conf_code]),
                 null_resp1_t_unconf = mean(d_block$resp1_t[d_block$stim_type == null_code & d_block$resp2 == unconf_code]),
                 null_resp1_t_no_conf = mean(d_block$resp1_t[d_block$stim_type == null_code& d_block$resp1 == no_code & d_block$resp2 == conf_code]),
                 null_resp1_t_no_unconf = mean(d_block$resp1_t[d_block$stim_type == null_code& d_block$resp1 == no_code & d_block$resp2 == unconf_code]),
                 
                 null_resp2_t = mean(d_block$resp2_t[d_block$stim_type == null_code]),
                 null_resp2_t_no = mean(d_block$resp2_t[d_block$stim_type == null_code & d_block$resp1 == no_code]),
                 null_resp2_t_conf = mean(d_block$resp2_t[d_block$stim_type == null_code & d_block$resp2 == conf_code]),
                 null_resp2_t_unconf = mean(d_block$resp2_t[d_block$stim_type == null_code & d_block$resp2 == unconf_code]),
                 null_resp2_t_no_conf = mean(d_block$resp2_t[d_block$stim_type == null_code& d_block$resp1 == no_code & d_block$resp2 == conf_code]),
                 null_resp2_t_no_unconf = mean(d_block$resp2_t[d_block$stim_type == null_code& d_block$resp1 == no_code & d_block$resp2 == unconf_code]),
                 
                 null_n = sum(d_block$stim_type == null_code), #, na.rm = TRUE
                 null_n_no = sum(d_block$stim_type == null_code & d_block$resp1 == no_code),
                 null_n_yes = sum(d_block$stim_type == null_code & d_block$resp1 == yes_code),
                 null_n_conf = sum(d_block$stim_type == null_code & d_block$resp2 == conf_code),
                 null_n_unconf = sum(d_block$stim_type == null_code & d_block$resp2 == unconf_code),
                 null_n_no_conf = sum(d_block$stim_type == null_code & d_block$resp1 == no_code & d_block$resp2 == conf_code),
                 null_n_no_unconf = sum(d_block$stim_type == null_code & d_block$resp1 == no_code & d_block$resp2 == unconf_code),
                 null_n_yes_conf = sum(d_block$stim_type == null_code & d_block$resp1 == yes_code & d_block$resp2 == conf_code),
                 null_n_yes_unconf = sum(d_block$stim_type == null_code & d_block$resp1 == yes_code & d_block$resp2 == unconf_code),
                 
                 near_trial_t = mean(d_block$trial_t[d_block$stim_type == near_code], na.rm = TRUE),
                 near_trial_t_no = mean(d_block$trial_t[d_block$stim_type == near_code & d_block$resp1 == no_code], na.rm = TRUE),
                 near_trial_t_yes = mean(d_block$trial_t[d_block$stim_type == near_code & d_block$resp1 == yes_code], na.rm = TRUE),
                 near_intensity = mean(d_block$intensity[d_block$stim_type == near_code]),
                 near_resp1 = mean(d_block$resp1[d_block$stim_type == near_code]),
                 near_resp1_unconf = mean(d_block$resp1[d_block$stim_type == near_code & d_block$resp2 == unconf_code]),
                 near_resp1_conf = mean(d_block$resp1[d_block$stim_type == near_code & d_block$resp2 == conf_code]),
                 near_resp2 = mean(d_block$resp2[d_block$stim_type == near_code]),
                 near_resp2_no = mean(d_block$resp2[d_block$stim_type == near_code & d_block$resp1 == no_code]),
                 near_resp2_yes = mean(d_block$resp2[d_block$stim_type == near_code & d_block$resp1 == yes_code]),
                 near_n = sum(d_block$stim_type == near_code),
                 near_n_yes = sum(d_block$stim_type == near_code & d_block$resp1 == yes_code),
                 near_n_no = sum(d_block$stim_type == near_code & d_block$resp1 == no_code),
                 near_n_conf = sum(d_block$stim_type == near_code & d_block$resp2 == conf_code),
                 near_n_unconf = sum(d_block$stim_type == near_code & d_block$resp2 == unconf_code),
                 near_n_yes_conf = sum(d_block$stim_type == near_code & d_block$resp1 == yes_code & d_block$resp2 == conf_code),
                 near_n_yes_unconf = sum(d_block$stim_type == near_code & d_block$resp1 == yes_code & d_block$resp2 == unconf_code),
                 near_n_no_conf = sum(d_block$stim_type == near_code & d_block$resp1 == no_code & d_block$resp2 == conf_code),
                 near_n_no_unconf = sum(d_block$stim_type == near_code & d_block$resp1 == no_code & d_block$resp2 == unconf_code),
                 
                 near_trial_median_yes = median(d_block$trial[d_block$stim_type == near_code & d_block$resp1 == yes_code]),
                 
                 near_resp1_t = mean(d_block$resp1_t[d_block$stim_type == near_code]),
                 near_resp1_t_no = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp1 == no_code]),
                 near_resp1_t_no_unconf = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp1 == no_code & d_block$resp2 == unconf_code]),
                 near_resp1_t_no_conf = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp1 == no_code & d_block$resp2 == conf_code]),
                 near_resp1_t_yes = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp1 == yes_code]),
                 near_resp1_t_yes_unconf = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp1 == yes_code & d_block$resp2 == unconf_code]),
                 near_resp1_t_yes_conf = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp1 == yes_code & d_block$resp2 == conf_code]),
                 near_resp1_t_conf = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp2 == conf_code]),
                 near_resp1_t_unconf = mean(d_block$resp1_t[d_block$stim_type == near_code & d_block$resp2 == unconf_code]),
                 
                 near_resp2_t = mean(d_block$resp2_t[d_block$stim_type == near_code]),
                 near_resp2_t_no = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp1 == no_code]),
                 near_resp2_t_no_unconf = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp1 == no_code & d_block$resp2 == unconf_code]),
                 near_resp2_t_no_conf = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp1 == no_code & d_block$resp2 == conf_code]),
                 near_resp2_t_yes = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp1 == yes_code]),
                 near_resp2_t_yes_unconf = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp1 == yes_code & d_block$resp2 == unconf_code]),
                 near_resp2_t_yes_conf = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp1 == yes_code & d_block$resp2 == conf_code]),
                 near_resp2_t_conf = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp2 == conf_code]),
                 near_resp2_t_unconf = mean(d_block$resp2_t[d_block$stim_type == near_code & d_block$resp2 == unconf_code]))
    
    
    # Update data frame
    db[nrow(db)+1,] <- results
    
    # If first row, update column names of data frame
    if (nrow(db)==1) { colnames(db) <- names(results)} 
  }
}


# BLOCK FILTER ------------------------------------------------------------

# Number of trials filter
db$n_filter <- as.integer(db$null_n >= null_n_min & db$near_n >= near_n_min)

# Detection filter
db$null_resp1_filter <- as.integer(db$null_resp1 <= null_resp1_max)
db$near_resp1_filter <- as.integer(db$near_resp1 >= near_resp1_min & db$near_resp1 <= near_resp1_max)

db$resp1_filter <- as.integer(db$null_resp1_filter & db$near_resp1_filter)

# Median near hit filter
db$near_trial_median_yes_filter <- as.integer(db$near_trial_median_yes > near_trial_median_yes_min & db$near_trial_median_yes <= near_trial_median_yes_max)

# More hit than false alarms?
db$HR_larger_FAR_filter <- as.integer((db$near_resp1-db$null_resp1)>margin_hit_FA)

# Combine filter
#db$block_filter <- as.integer(db$n_filter & db$resp1_filter & db$near_trial_median_yes_filter)
db$block_filter <- as.integer(db$n_filter & db$resp1_filter & db$HR_larger_FAR_filter)


return(db)

}
