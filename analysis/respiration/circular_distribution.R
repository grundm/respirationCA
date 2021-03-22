###########################################################################
####                       Circular Distribution                       ####
###########################################################################

# Author:         Martin Grund, Marc Pabst
# Last update:    February 9, 2020


# Settings ---------------------------------------------------------------

rm(list=ls())

code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'


# Load data ---------------------------------------------------------------

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210208.Rds", sep="/"))


# Select valid data -------------------------------------------------------

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_filter == 1)


# Plot circular distribution within respiratory cycle ---------------------

source(paste(code_dir, '/ecg/plot_test_R.R', sep = ''))

output_pdf <- '~/ownCloud/promotion/experiment/respirationCA/fig/resp_circ_20210209.pdf'

pdf(output_pdf, width=6, height=6, useDingbats=FALSE)

angle_hit <- plot_test_R(data_stims_valid,1,1,c(0,1),'Hit')

angle_miss <- plot_test_R(data_stims_valid,1,0,c(0,1),'Miss')

angle_CR <- plot_test_R(data_stims_valid,0,0,c(0,1),'CR')

angle_FA <- plot_test_R(data_stims_valid,0,1,c(0,1),'FA')

angle_hit_conf <- plot_test_R(data_stims_valid,1,1,1,'Hit conf')

angle_hit_unconf <- plot_test_R(data_stims_valid,1,1,0,'Hit unconf')

angle_miss_conf <- plot_test_R(data_stims_valid,1,0,1,'Miss conf')

angle_miss_unconf <- plot_test_R(data_stims_valid,1,0,0,'Miss unconf')

angle_CR_conf <- plot_test_R(data_stims_valid,0,0,1,'CR conf')

angle_CR_unconf <- plot_test_R(data_stims_valid,0,0,0,'CR unconf')

dev.off()

# Post-processing:
system(paste('mkdir ', dirname(output_pdf), '/resp_circ', sep = ''))
system(paste('convert -density 300 ', output_pdf,' -colorspace RGB ', dirname(output_pdf), '/resp_circ/resp_circ-%d.png', sep = ''))


# Moore's test ------------------------------------------------------------

source(paste(code_dir, '/respiration/circular_tests.R', sep = ''))

AngleStats <- MooreRStats(conversion.circular(angle_hit$angle, type="angles", units="radians", rotation="clock", zero=0),
                          conversion.circular(angle_miss$angle, type="angles", units="radians", rotation="clock", zero=0))

cosphi <- AngleStats[[1]][!is.na(AngleStats[[1]])]
sinphi <- AngleStats[[2]][!is.na(AngleStats[[1]])]
Ranks <- AngleStats[[3]][!is.na(AngleStats[[1]])]
NR <- 10000

MooreRTestRand(cosphi, sinphi, Ranks, NR)


# Angular variance --------------------------------------------------------

t.test(angle_hit$var_angle,angle_miss$var_angle,paired = T)
t.test(angle_hit$var_angle,angle_CR$var_angle,paired = T)

t.test(angle_miss$var_angle,angle_CR$var_angle,paired = T)


# Entrainment? ------------------------------------------------------------

# Binning trials
intervals <- seq(0,150,20)

for (i in 1:(length(intervals)-1)) {
  
  interval_start_tmp <- intervals[i]
  interval_end_tmp <- intervals[i+1]
  
  data_stims_valid_tmp <- subset(data_stims_valid, block == 1 & trial > interval_start_tmp & trial <= interval_end_tmp)
  
  angle_tmp <- plot_test_R(data_stims_valid_tmp,1,c(0,1),c(0,1),paste('All: ', interval_start_tmp, '-', interval_end_tmp, ' trial', sep = ''))
}


# Detection rate across respiratory cycle ---------------------------------

# Bin and average hit rates -----------------------------------------------

# Distance of stimulus onset to the previous R peak
data_stims$dist=NaN

# 4 x 1-s interval 0-4 s
intervals <- seq(0,4,1)

for (i in 1:(length(intervals)-1)) {
  print(i)
  print(intervals[i])
  print(intervals[i+1])
  data_stims$dist[data_stims$diff2inhale_onset>=intervals[i] & data_stims$diff2inhale_onset<intervals[i+1]] <- i
}

# Select only the near-threshold trials
data_select = subset(data, stim_type==1 & !is.nan(dist) & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1)

# Average detection for each participant and bin relative to stimulus onset (see above)
HR=aggregate(resp1 ~ ID*dist, data_select, FUN=mean)

HR$ID <- as.factor(HR$ID)
HR$dist <- as.factor(HR$dist)

library("afex")

fit_all <- aov_ez("ID","resp1", HR, within="dist")
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values
