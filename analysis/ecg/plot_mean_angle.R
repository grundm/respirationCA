# respirationCA - ECG plot mean angle for hits and misses -----

# Author:         Martin Grund
# Last update:    July 29, 2021


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

# Trial type variable -----------------------------------------------------

data$trial_type <- NA

data$trial_type[data$stim_type==0 & data$resp1==0] <- 'CR'
data$trial_type[data$stim_type==0 & data$resp1==1] <- 'FA'
data$trial_type[data$stim_type==1 & data$resp1==0] <- 'near_miss'
data$trial_type[data$stim_type==1 & data$resp1==1] <- 'near_hit'

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

# Exclude false alarms
#data_valid <- subset(data_valid, trial_type != 'FA')

par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

res<-hist(data_valid$stim_degree[data_valid$stim_type == 1],
     main = 'Stimulus onsets relative to cardiac cycle\n across all near-threshold trails and participants',
     xlab = 'Degree within cardiac cycle (0 = R-peak)',
     ylab = 'Cumulative number of trials',
     xaxt = 'n',
     ylim = c(0,900))

axis(side=1, at=seq(0,360,20), labels=seq(0,360,20))

rtest <- rayleigh.test(circular(data_valid$stim_degree[data_valid$stim_type == 1], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0))

rtest$statistic
rtest$p.value

angle_near <- plot_test_R(data_valid,1,c(0,1),c(0,1),'All near-threshold trials')

min(p.adjust(angle_near$r_pval, 'fdr'))
#angle_near$r_pval_fdr <- p.adjust(angle_near$r_pval, 'fdr')

# Check all participants --------------------------------------------------

source(paste(code_dir, '/ecg/plot_test_R.R', sep = ''))

# Exclude false alarms
data_valid <- subset(data_valid, trial_type != 'FA')

angle_all_trials <- plot_test_R(data_valid,c(0,1),c(0,1),c(0,1),'All trials')

angle_all_trials$r_pval_fdr <- p.adjust(angle_all_trials$r_pval, 'fdr')

par(mfrow=c(7,6))

for (i in unique(data_valid$ID)) {
  
  #plot_test_R(subset(data_valid, ID == i),c(0,1),c(0,1),c(0,1),paste('ID', i, ' - all trials', sep = ''))
  
  angle <- circular(data_valid$stim_degree[data_valid$ID==i], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=pi/2)
  #angle <- circular(data_valid$stim_degree, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=pi/2)
  
  # Get FDR-corrected Rayleigh test p-value
  r_pval_fdr_tmp <- angle_all_trials$r_pval_fdr[angle_all_trials$ID_list == i]
  
  plot_col <- ifelse(r_pval_fdr_tmp < 0.05, 'orange', 'grey25')
  
  # default par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  par(mar = c(0.5, 0, 1, 1))
  plot(mean(angle,na.rm=T),
       stack=FALSE,
       axes=FALSE,
       bins = 180,
       col = "gray35")#,
       #cex = 1.6)

  title(i)#, line = 4, cex.main = 2)
  
  # Arrow 
  #arrows.circular(angle, col = "grey80", lwd = 2)
  #arrows.circular(angle, y = mean_length, col = "grey80", lwd = 2, code = 0)
  arrows.circular(mean(angle,na.rm=T), y=rho.circular(angle, na.rm=T), lwd=1, col = plot_col)
  circ.dens = density(angle[!is.na(angle)], bw=20)
  #lines(circ.dens, col="grey25", lwd = 3, xpd=TRUE, shrink = 2) 
  lines(circ.dens, col = plot_col, lwd = 2, xpd=TRUE) 
  
  # Rayleigh test for uniformity
  r_test <- rayleigh.test(angle)
  #print(r_test)
  
  text(0,-.25,paste('R = ', round(r_test$statistic,2), sep = ''),cex = 1.1)
  text(0,-.65,paste('p = ', signif(r_pval_fdr_tmp,2), sep = ''),cex = 1.1)
}

# Plot each stimulus-response condition -----------------------------------

source(paste(code_dir, '/ecg/plot_test_R.R', sep = ''))

output_pdf <- '~/ownCloud/promotion/experiment/respirationCA/fig/cardiac_circ_20210721_v2.pdf'

pdf(output_pdf, width=6, height=6, useDingbats=FALSE)

# All trials
plot_test_R(data_valid,c(0,1),c(0,1),c(0,1),'All trials')

#angle_conf <- plot_test_R(data_valid,c(0,1),c(0,1),1,'All confident trials')
#angle_conf <- plot_test_R(data_valid,c(0,1),c(0,1),0,'All unconfident trials')

# All near-threshold trials
angle_near <- plot_test_R(data_valid,1,c(0,1),c(0,1),'All near-threshold trials')

# All confident trials w/o false alarms
angle_conf <- plot_test_R(subset(data_valid, trial_type != 'FA'),c(0,1),c(0,1),1,'All confident trials')
angle_unconf <- plot_test_R(subset(data_valid, trial_type != 'FA'),c(0,1),c(0,1),0,'All unconfident trials')

# Near hit (R = 0.34, p = 0.007)
angle_hit <- plot_test_R(data_valid,1,1,c(0,1),'Hit')

# Near miss (R = 0.20, p = 0.18)
angle_miss <- plot_test_R(data_valid,1,0,c(0,1),'Miss')

# CR (R = 0.05, p = 0.91)
angle_CR <- plot_test_R(data,0,0,c(0,1),'Correct rejection')

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

# Post-processing:
output_dir <- paste(dirname(output_pdf), '/', tools::file_path_sans_ext(basename(output_pdf)), sep = '')
system(paste('mkdir ', output_dir, sep = ''))
system(paste('convert -density 300 ', output_pdf,' -colorspace RGB ', output_dir, '/cardiac_circ-%d.png', sep = ''))


# Moore's test ------------------------------------------------------------

source(paste(code_dir, '/respiration/circular_tests.R', sep = ''))

AngleStats <- MooreRStats(conversion.circular(angle_miss$angle, type="angles", units="radians", rotation="clock", zero=0),
                          conversion.circular(angle_CR$angle, type="angles", units="radians", rotation="clock", zero=0))

cosphi <- AngleStats[[1]][!is.na(AngleStats[[1]])]
sinphi <- AngleStats[[2]][!is.na(AngleStats[[1]])]
Ranks <- AngleStats[[3]][!is.na(AngleStats[[1]])]
NR <- 10000

MooreRTestRand(cosphi, sinphi, Ranks, NR)


# Confidence GLMM ---------------------------------------------------------

# Analysis of the effect of the physiological measures on confidence,
# while controlling for task performance, i.e.,
# confidence ~ respiratory phase + behav_resp, 
# and confidence ~ cardiac phase + behav_resp 
# (where behav_resp is a binary regressor "hit vs miss").

# Kees Mulder (https://stats.stackexchange.com/users/55410/kees-mulder), Hierarchical circular-circular regression in R, URL (version: 2021-03-18): https://stats.stackexchange.com/q/514198

# Transform degrees to radians
data_valid$stim_degree_rad <- data_valid$stim_degree * pi / 180
#data_valid$stim_degree_rad2 <- circular(data_valid$stim_degree_rad, type="angles", units="radians", rotation="clock", zero=0)

# # null model (no effect of cycle)
# model_null = glmer(resp1 ~ 1 + (1|ID),
#                    family = 'binomial',
#                    data = data_valid)
# 
# # additive model (effect of cycle but no interaction)
# model_additive = glmer(resp1 ~ 1 + cos(stim_degree_rad) + sin(stim_degree_rad) + (1|ID),
#                        family = 'binomial',
#                        data = data_valid)

# null model (no effect of cycle)
model_null = glmer(resp2 ~ 1 + (1|ID),
                   family = 'binomial',
                   data = data_valid)

# additive model (effect of cycle but no interaction)
model_additive = glmer(resp2 ~ 1 + cos(stim_degree_rad) + sin(stim_degree_rad) + (1|ID),
                       family = 'binomial',
                       data = data_valid)

anova(model_null,model_additive)


# null model (no effect of cycle)
model_null = glmer(resp2 ~ resp1 + (1|ID),
                   family = 'binomial',
                   data = subset(data_valid, stim_type == 1))

# additive model (effect of cycle but no interaction)
model_additive = glmer(resp2 ~ resp1 + cos(stim_degree_rad) + sin(stim_degree_rad) + (1|ID),
                       family = 'binomial',
                       data = subset(data_valid, stim_type == 1))

anova(model_null,model_additive)

model_additive = glmer(resp2 ~ resp1 + diff2peak + (1|ID),
                       family = 'binomial',
                       data = data_valid)
