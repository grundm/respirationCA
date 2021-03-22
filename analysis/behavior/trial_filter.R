function(d) {
# Preprocessing on trial level

# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      November 12, 2020

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

cue_t_min <- 0.61
cue_t_max <- 0.62

# TRIAL DURATION ----------------------------------------------------------

d$trial_t <- ifelse(d$trial < max(d$trial), tail(d$onset_fix,-1) - head(d$onset_fix,-1), d$onset_resp2-d$onset_fix+resp2_max)

d$trial_t2 <- ifelse(d$trial < max(d$trial), tail(d$onset_cue,-1) - head(d$onset_cue,-1), NA)

d$ISI[2:nrow(d)] <- d$stim_onset[2:nrow(d)] - d$stim_onset[1:(nrow(d)-1)]

d$ISI[d$trial == 1] <- NA;


# FILTER ----------------------------------------------------

# Response time filter
d$resp_t_filter <- as.integer(d$resp1_t > resp1_min & d$resp1_t < resp1_max & d$resp2_t > resp2_min & d$resp2_t < resp2_max)

# Response button filter
d$resp_btn_filter <- as.integer(match(d$resp1_btn, resp1_btn, 0) & match(d$resp2_btn, resp2_btn, 0))

# Cue presentation time filter
d$cue_t_filter <- as.integer( (d$onset_resp1 - d$onset_cue) > cue_t_min & (d$onset_resp1 - d$onset_cue) < cue_t_max)

# Combine response time & button filter
#d$resp_filter <- as.integer(d$resp_t_filter & d$resp_btn_filter & d$cue_t_filter)
d$resp_filter <- as.integer(d$resp_t_filter & d$resp_btn_filter)

return(d)

}
