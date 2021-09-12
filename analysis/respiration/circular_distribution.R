###########################################################################
####                       Circular Distribution                       ####
###########################################################################

# Author:         Martin Grund, Marc Pabst
# Last update:    August 6, 2021


# Settings ---------------------------------------------------------------

rm(list=ls())

code_dir <- '~/ownCloud/promotion/experiment/respirationCA/code/respirationca'
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'


# Load data ---------------------------------------------------------------

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210622.Rds", sep="/")) # Locked to inspiration and expiration

# Select valid data -------------------------------------------------------

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_exhale_filter == 1)
data_stims_valid$stim_degree <- data_stims_valid$stim_degree_exhale

#data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_inhale_filter == 1)
#data_stims_valid$stim_degree <- data_stims_valid$stim_degree_inhale

# Exclude false alarms ----------------------------------------------------

data_stims_valid$FA_filter <- ifelse(data_stims_valid$stim_type == 0 & data_stims_valid$resp1 == 1, 1, 0)
sum(data_stims_valid$FA_filter)
data_stims_valid <- subset(data_stims_valid, FA_filter == 0)

# Check all participants --------------------------------------------------

par(mfrow = c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

bin_res <- hist(data_stims_valid$stim_degree,
                #main = 'Stimulus onsets relative to respiration cycle\n across all trails and participants',
                main = 'Trials across respiration cycle',
                #xlab = 'Degree within respiration cycle (0 = expiration onset)',
                xlab = 'Respiration phase',
                ylab = 'Number of trials',
                xaxt = 'n',
                #ylim = c(0,2000),
                col = rgb(0,0.447,0.741,0.7),
                breaks = seq(0,360,20))#,
                #border = F)

axis(side=1, at=seq(0,360,20), labels=paste(seq(0,360,20),'°',sep=''))

label_list <- paste(as.character(head(bin_res$breaks,-1)),'\n-',as.character(bin_res$breaks[-1]),'°',sep = '')
names(label_list) <- 1:(length(bin_res$breaks)-1)

axis(side=1, at=head(bin_res$breaks,-1)+(0.5*bin_res$breaks[2]), labels=label_list, lwd = 0, cex.axis = 0.7)

rtest <- rayleigh.test(circular(data_stims_valid$stim_degree, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0))

rtest$statistic
rtest$p.value

breaks_tmp <- bin_res$breaks[2:length(bin_res$breaks)]
breaks_tmp[breaks_tmp > 320 | breaks_tmp <= 100]

sum(bin_res$counts[breaks_tmp > 320 | breaks_tmp <= 100])/sum(bin_res$counts)

sum(bin_res$counts[breaks_tmp > 100 & breaks_tmp <= 320])/sum(bin_res$counts)

# Check all participants --------------------------------------------------

angle_all_trials <- plot_test_R(data_stims_valid,c(0,1),c(0,1),c(0,1),'All trials')

angle_all_trials$r_pval_fdr <- p.adjust(angle_all_trials$r_pval, 'fdr')

par(mfrow=c(7,6))

for (i in unique(data_stims_valid$ID)) {
  
  # Get angles of each participants
  angle <- circular(data_stims_valid$stim_degree[data_stims_valid$ID==i], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=pi/2)
  #angle <- circular(data_valid$stim_degree, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=pi/2)
  
  # Get FDR-corrected Rayleigh test p-value
  r_pval_fdr_tmp <- angle_all_trials$r_pval_fdr[angle_all_trials$ID_list == i]
  
  plot_col <- ifelse(r_pval_fdr_tmp < 0.05, 'orange', 'grey25')
  
  # default par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  par(mar = c(1, 0, 1, 1))
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
  lines(circ.dens, col=plot_col, lwd = 2, xpd=TRUE) 
  
  # Rayleigh test for uniformity
  r_test <- rayleigh.test(angle)
  
  text(0,-.25,paste('R = ', round(r_test$statistic,2), sep = ''),cex = 1.1)
  text(0,-.65,paste('p = ', signif(r_pval_fdr_tmp,2), sep = ''),cex = 1.1)
}

# Plot circular distribution within respiratory cycle ---------------------

source(paste(code_dir, '/ecg/plot_test_R.R', sep = ''))

output_pdf <- '~/ownCloud/promotion/experiment/respirationCA/fig/resp_circ_20210716.pdf'

pdf(output_pdf, width=6, height=6, useDingbats=FALSE)

angle_all_trials <- plot_test_R(data_stims_valid,c(0,1),c(0,1),c(0,1),'All trials')

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
output_dir <- paste(dirname(output_pdf), '/', tools::file_path_sans_ext(basename(output_pdf)), sep = '')
system(paste('mkdir ', output_dir, sep = ''))
system(paste('convert -density 300 ', output_pdf,' -colorspace RGB ', output_dir, '/resp_circ-%d.png', sep = ''))


# Does variance correlate with hit rate -----------------------------------

ggscatter(angle_all_trials, x = "var_angle", y = "HR",
          add = 'reg.line',
          cor.coef = TRUE,
          cor.coef.coord = c(0.40,0.80)) +
  #stat_regline_equation(label.y = 0.99) +
  stat_regline_equation(label.x = 0.40, label.y = 0.75, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +
  labs(x = "Circular variance of all trials", y = "Near-threshold hit rate")


ggscatter(angle_CR, x = "var_angle", y = "HR",
          add = 'reg.line',
          cor.coef = TRUE) +
  #stat_regline_equation(label.y = 0.99) +
  stat_regline_equation(label.y = 0.75, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +
  labs(x = "Circular variance for correct rejections", y = "Near-threshold hit rate")


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


# Inspiration onset in respiration cycle (relative exhale onset) ----------

data_stims_valid$inhale_onset_degree_exhale <- (data_stims_valid$exhale_duration_exhale/data_stims_valid$resp_cycle_t_exhale)*360

# data_stims_valid$stim_degree <- data_stims_valid$stim_degree_exhale
data_stims_valid$stim_degree <- data_stims_valid$inhale_onset_degree_exhale

source(paste(code_dir, '/ecg/plot_test_R.R', sep = ''))

angle_inhale_onset_exhale <- plot_test_R(data_stims_valid,c(0,1),c(0,1),c(0,1),'Inhale onset')

#angle_tmp <- circular(data_stims_valid$inhale_onset_degree_exhale, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=pi/2)
#plot(angle_tmp)
#arrows.circular(mean(angle_tmp,na.rm=T), y=rho.circular(angle_tmp, na.rm=T), lwd=1, col = "grey25")

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
# 4 x 90° interval 0-360°
intervals <- seq(0,360,90)
intervals <- seq(0,360,45)

for (i in 1:(length(intervals)-1)) {
  print(i)
  print(intervals[i])
  print(intervals[i+1])
  #data_stims$dist[data_stims$diff2inhale_onset>=intervals[i] & data_stims$diff2inhale_onset<intervals[i+1]] <- i
  #data_stims$dist[data_stims$diff2exhale_onset>=intervals[i] & data_stims$diff2exhale_onset<intervals[i+1]] <- i
  data_stims$dist[data_stims$stim_degree_exhale>=intervals[i] & data_stims$stim_degree_exhale<intervals[i+1]] <- i
}

# Select only the near-threshold trials
data_select = subset(data_stims, stim_type==1 & !is.nan(dist) & resp_filter == 1 & HR_larger_FAR_filter == 1 & resp_cycle_exhale_filter == 1)

# Merge 1st and 4th interval: data_select$dist[data_select$dist==4] = 1
# Merge 2nd and 3rd interval: data_select$dist[data_select$dist==3] = 2

# Average detection for each participant and bin relative to stimulus onset (see above)
HR = aggregate(resp1 ~ ID*dist, data_select, FUN=mean)

HR$ID <- as.factor(HR$ID)
HR$dist <- as.factor(HR$dist)

library("afex")

fit_all <- aov_ez("ID","resp1", HR, within="dist")
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values

boxplot(resp1 ~ dist, HR, 
        ylab = 'Hit rate',
        #names = c('270-90°','90-270°'))
        names = c('0-90°','90-180°','180-270°','270-360°'))
title('Near-threshold trials binned in respiration cycle\n (start expiration)')
aggregate(resp1~dist,HR,FUN=mean)

library(rstatix)
library(dplyr)

stat_test = subset(HR, ID != 6) %>%
  t_test(resp1 ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE)

stat_test_HR = HR %>%
  t_test(resp1 ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.7)

label_list <- paste(as.character(head(intervals,-1)*1000),as.character(intervals[-1]*1000),sep = "\n-")
label_list <- paste(as.character(head(intervals,-1)),'-',as.character(intervals[-1]),'°',sep = '')
names(label_list) <- 1:(length(intervals)-1)

# Plot hit-rates for bins
HR %>%
  ggplot(aes(y=resp1, x=dist)) + 
  geom_boxplot(fill = rgb(251,128,114, maxColorValue=255)) +
  scale_x_discrete(labels=label_list) +
  ggtitle("Hit rate x respiration phase") +
  ylab("Hit rate") + 
  xlab("Respiration phase") +
  theme_cowplot() +
  theme(legend.title = element_blank()) +
  #scale_fill_manual(rgb(251,128,114, maxColorValue=255)) +
  stat_pvalue_manual(stat_test_HR, label = "p.adj", tip.length = 0.01, hide.ns=T)
