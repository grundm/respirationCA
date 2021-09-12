# respirationCA - Plot dprime,  criterion,  and m ratios ------------------


# Author:         Martin Grund
# Last update:    August 5, 2021


# Load data ---------------------------------------------------------------

rm(list = ls())

data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data/Mratio_new'

library(R.matlab)

fit01 <- readMat(paste(data_dir, '/fit_first.mat', sep = ''))
fit02 <- readMat(paste(data_dir, '/fit_second.mat', sep = ''))
fit03 <- readMat(paste(data_dir, '/fit_third.mat', sep = ''))
fit04 <- readMat(paste(data_dir, '/fit_fourth.mat', sep = ''))

library(dplyr)
library(rstatix)
library(ggpubr)
library(cowplot)

# dprime ------------------------------------------------------------------

d <- data.frame()

d <- cbind(matrix(1,41,1), t(fit01$fit.first[[3]]))
d <- rbind(d, cbind(matrix(2,41,1), t(fit02$fit.second[[3]])))
d <- rbind(d, cbind(matrix(3,41,1), t(fit03$fit.third[[3]])))
d <- rbind(d, cbind(matrix(4,41,1), t(fit04$fit.fourth[[3]])))

d <- data.frame(d)
colnames(d) <- c('int', 'dprime')
d$int <- as.factor(d$int)

stat_test_d = d %>%
  t_test(dprime ~ int, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

d %>%
  ggplot(aes(y=dprime, x=int)) + 
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_brewer(palette="Dark2")  +
  stat_pvalue_manual( stat_test_d, label = "p", tip.length = 0.01, hide.ns=T)+
  theme(legend.position="none")




# M Ratio -----------------------------------------------------------------

mr <- data.frame()
# 11 - M-ratio S1; 12 - M-ratio S2
mr <- cbind(c(1:40, 1001), matrix(1,41,1), t(fit01$fit.first[[11]]), t(fit01$fit.first[[12]]))
mr <- rbind(mr, cbind(c(1:40, 1001), matrix(2,41,1), t(fit02$fit.second[[11]]), t(fit02$fit.second[[12]])))
mr <- rbind(mr, cbind(c(1:40, 1001), matrix(3,41,1), t(fit03$fit.third[[11]]), t(fit03$fit.third[[12]])))
mr <- rbind(mr, cbind(c(1:40, 1001), matrix(4,41,1), t(fit04$fit.fourth[[11]]), t(fit04$fit.fourth[[12]])))

mr <- data.frame(mr)
colnames(mr) <- c('ID', 'dist', 'no', 'yes')
mr$dist <- as.factor(mr$dist)

mr_S1_S2 <- pivot_longer(mr, cols = c(yes,no), names_to = 'S1_S2', values_to = 'mratio')

# Plot M-ratios for no and yes-response
mr_S1_S2 %>%
  ggplot(aes(y=mratio, x=dist, fill=S1_S2)) + 
  geom_boxplot() +
  scale_x_discrete(labels=c("1" = "0-200", "2" = "200-400", "3" = "400-600", "4" = "600-800")) +
  ggtitle("M-ratio x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("M-ratio") +
  theme_cowplot() +
  theme(legend.title = element_blank()) +
  #scale_fill_manual(values = c('#e78ac3', '#66c2a5'))
  scale_fill_manual(values = c('#e78ac3', '#66c2a5'))
  #scale_fill_manual(values = c('#66c2a5', '#fc8d62'))
  #scale_fill_brewer(palette="Dark2")  +
  #stat_pvalue_manual( stat_test_inhale, label = "p", tip.length = 0.01, hide.ns=T)+
  #theme(legend.position="none")

library("afex")

fit_all <- aov_ez("ID","mratio", mr_S1_S2, within=c("S1_S2", "dist"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values

aggregate(mratio ~ dist*S1_S2, mr_S1_S2, FUN='mean')

# Plot M-ratios separately for yes and no ---------------------------------

stat_test_mr_yes <- subset(mr_S1_S2, S1_S2 == 'yes') %>%
  t_test(mratio ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

# Plot M-ratios for yes-response
subset(mr_S1_S2, S1_S2 == 'yes') %>%
  ggplot(aes(y=mratio, x=dist)) + 
  geom_boxplot() +
  scale_x_discrete(labels=c("1" = "0-200", "2" = "200-400", "3" = "400-600", "4" = "600-800")) +
  ggtitle("M-ratio (yes) x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("M-ratio") +
  theme_cowplot() +
  stat_pvalue_manual(stat_test_mr_yes, label = "p.adj", tip.length = 0.01, hide.ns=T)
#scale_fill_brewer(palette="Dark2")  +
#stat_pvalue_manual(stat_test_mr_yes, label = "p.adj", tip.length = 0.01, hide.ns=T)+
#theme(legend.position="none")


stat_test_mr_no <- subset(mr_S1_S2, S1_S2 == 'no') %>%
  t_test(mratio ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.3)

# Plot M-ratios for no-response
subset(mr_S1_S2, S1_S2 == 'no') %>%
  ggplot(aes(y=mratio, x=dist)) + 
  geom_boxplot() +
  scale_x_discrete(labels=c("1" = "0-200", "2" = "200-400", "3" = "400-600", "4" = "600-800")) +
  ggtitle("M-ratio (no) x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("M-ratio") +
  theme_cowplot() +
  stat_pvalue_manual(stat_test_mr_no, label = "p.adj", tip.length = 0.01, hide.ns=T)

p.adjust(c(stat_test_mr_no$p, stat_test_mr_yes$p),'fdr')<0.05

# Plot false alarms across cardiac cycle ----------------------------------

# Create labels based on intervals
label_list <- paste(as.character(head(intervals,-1)*1000),as.character(intervals[-1]*1000),sep = "\n-")
names(label_list) <- 1:(length(intervals)-1)

null_near_n$dist <- as.factor(null_near_n$dist)

stat_test_FA <- subset(null_near_n, stim_type != 0) %>%
  t_test(yes_rate ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>% 
  add_xy_position(step.increase=0.5)

# Plot number of false alarms
subset(null_near_n, stim_type != 0) %>%
  ggplot(aes(y=yes_rate, x=dist)) + 
  geom_boxplot() +
  scale_x_discrete(labels=label_list) +
  ggtitle("Hit rate x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("Hit rate") +
  theme_cowplot() +
  stat_pvalue_manual(stat_test_FA, label = "p.adj", tip.length = 0.01, hide.ns=T)


# Calculate dprime and c --------------------------------------------------

data_select = subset(data, stim_type==1 & !is.nan(dist) & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1 & ID != 6)
HR = aggregate(resp1 ~ ID*dist, data_select, FUN=mean)
near_yes = aggregate(resp1 ~ ID*dist, data_select, FUN=sum)


data_select = subset(data, stim_type==0 & !is.nan(dist) & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1 & ID != 6)
FAR = aggregate(resp1 ~ ID*dist, data_select, FUN=mean)



# d-prime and c -----------------------------------------------------------

# Select valid trials
data_select = subset(data, !is.nan(dist) & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1)

# Count null and near trials
null_near_n <- data_select %>% count(ID, dist, stim_type, name="n")

# Count only yes responses for null and near trials
yes_n <- subset(data_select, resp1 == 1) %>% count(ID, dist, stim_type, name="n")

# Add number of yes-responses to trial numbers
null_near_n$yes_n <- NA

for (i in c(1:nrow(null_near_n))) {
  if( length(yes_n$n[yes_n$ID == null_near_n$ID[i] & yes_n$dist == null_near_n$dist[i] & yes_n$stim_type == null_near_n$stim_type[i]]) != 0) {
    null_near_n$yes_n[i] <- yes_n$n[yes_n$ID == null_near_n$ID[i] & yes_n$dist == null_near_n$dist[i] & yes_n$stim_type == null_near_n$stim_type[i]]
  }
  else {
    null_near_n$yes_n[i] <- 0
  }
}

# Correct for zero false alarms (loglinear approach, Hautus (1995))
null_near_n$yes_rate <- (null_near_n$yes_n) / (null_near_n$n)
null_near_n$yes_rate_corr <- (null_near_n$yes_n + 0.5) / (null_near_n$n + 1.0)

# Long to wide format
library(tidyr)

FAR_HR <- null_near_n %>% pivot_wider(id_cols = c(ID,dist), names_from = stim_type, values_from = yes_rate_corr, names_prefix = 'yes_rate_corr_')

FAR_HR$dprime <- qnorm(FAR_HR$yes_rate_corr_1) - qnorm(FAR_HR$yes_rate_corr_0)

FAR_HR$c <- -(qnorm(FAR_HR$yes_rate_corr_1) + qnorm(FAR_HR$yes_rate_corr_0))/2

FAR_HR %>% t_test(dprime ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE)

FAR_HR$dist <- as.factor(FAR_HR$dist)

FAR_HR %>%
  ggplot(aes(y=dprime, x=dist)) + 
  geom_boxplot()

boxplot(dprime ~ dist, FAR_HR)

FAR_HR %>% t_test(c ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE)
boxplot(c ~ dist, FAR_HR)
