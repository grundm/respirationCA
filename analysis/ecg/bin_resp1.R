# respirationCA - Bin hit rates relative to R peak ------------------------

# Author:         Martin Grund, Esra Al
# Last update:    February 17, 2022

# Trials are assigned to bins depending on the stimulus onset relative to
# the previous R-peak. Hit rates in near-threshold trials are averaged for 
# each bin and participant. The results are plotted as boxplots and tested 
# with an ANOVA and post-hoc t-tests.

# Settings ----------------------------------------------------------------

rm(list=ls())

data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

data_file <- paste(data_dir, '/ecgdata_all.csv', sep = '')


# Read data ---------------------------------------------------------------

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

# Filter for hit rate > false alarm rate
for (i in unique(data$ID)) {
  for (b in unique(data$block[data$ID == i])) {
    HR <- mean(data$resp1[data$ID == i & data$block == b & data$resp_filter == 1 & data$stim_type == 1])
    FAR <- mean(data$resp1[data$ID == i & data$block == b & data$resp_filter == 1 & data$stim_type == 0])
    
    data$HR_larger_FAR_filter[data$ID == i & data$block == b] <- as.integer((HR-FAR)>margin_hit_FA)
  }
}


# Bin and average hit rates -----------------------------------------------

# Distance of stimulus onset to the previous R peak
data$dist=NaN

# 4 x 200-ms interval 0-800 ms
intervals <- seq(0,0.8,0.2)

# 40 x 20-ms interval 0-800 ms
intervals <- seq(0,0.8,0.02)

intervals <- seq(0,0.8,0.04)

intervals <- seq(0,0.8,0.05)
intervals <- seq(0,0.7,0.05)

for (i in 1:(length(intervals)-1)) {
  print(i)
  print(intervals[i])
  print(intervals[i+1])
  data$dist[data$diff2peak>=intervals[i] & data$diff2peak<intervals[i+1]] <- i
}

# Select only the near-threshold trials
data_select = subset(data, stim_type==1 & !is.nan(dist) & data$resp_filter == 1 & data$HR_larger_FAR_filter == 1)

# Average detection for each participant and bin relative to stimulus onset (see above)
HR=aggregate(resp1 ~ ID*dist, data_select, FUN=mean)
HR=aggregate(resp1 ~ ID*dist*resp2, data_select, FUN=mean)

HR$ID <- as.factor(HR$ID)
HR$dist <- as.factor(HR$dist)
HR$resp2 <- as.factor(HR$resp2)

# Plot hit rates across bins ----------------------------------------------

library(ggplot2)

# Boxplot of the hit rates for each time bin 
ggplot(HR, aes(x=dist, y=resp1, fill=resp2)) + # col=resp2
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Detection rate x confidence x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("Detection rate") +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position=c(0.12,0.10)
  ) +
  scale_x_discrete(labels=c("1" = "0-200", "2" = "200-400", "3" = "400-600", "4" = "600-900")) +
  scale_fill_discrete(labels=c("unconfident", "confident")) +
  annotate(geom = "text", x = 2.7, y = 0.08, label = "p = 0.02", color = "black") + # conf: T2 vs T3
  annotate("segment", x = 2.2, xend = 3.2, y = 0.05, yend = 0.05, colour = "black") +
  annotate(geom = "text", x = 3.2, y = -0.01, label = "p = 0.02", color = "black") +  # conf: T2 vs T4
  annotate("segment", x = 2.2, xend = 4.2, y = -0.04, yend = -0.04, colour = "black") +
  annotate(geom = "text", x = 1, y = 1.06, label = "p = 0.02", color = "black") +  # T1: conf vs unconf
  annotate("segment", x = 0.8, xend = 1.2, y = 1.03, yend = 1.03, colour = "black") +
  annotate(geom = "text", x = 2, y = 1.06, label = "p = 0.02", color = "black") +  # T2: conf vs unconf
  annotate("segment", x = 1.8, xend = 2.2, y = 1.03, yend = 1.03, colour = "black")

#+  stat_summary(aes(y = resp2,group=2), fun.y=mean, colour="red", geom="line",group=2)

#annotate(geom = "text", x = 3, y = 1.06, label = "p = 0.04", color = "black") +  # T3: conf vs unconf
#  annotate("segment", x = 2.8, xend = 3.2, y = 1.03, yend = 1.03, colour = "black")

# T-tests between bins ----------------------------------------------------

pval <- numeric(16)

# Conf vs. unconf within each interval
pval[1] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 0], HR$resp1[HR$dist==1 & HR$resp2 == 1], paired=T)$p.value # p = 0.003
pval[2] <- t.test(HR$resp1[HR$dist==2 & HR$resp2 == 0], HR$resp1[HR$dist==2 & HR$resp2 == 1], paired=T)$p.value # p = 0.002
pval[3] <- t.test(HR$resp1[HR$dist==3 & HR$resp2 == 0], HR$resp1[HR$dist==3 & HR$resp2 == 1], paired=T)$p.value # p = 0.04
pval[4] <- t.test(HR$resp1[HR$dist==4 & HR$resp2 == 0], HR$resp1[HR$dist==4 & HR$resp2 == 1], paired=T)$p.value # p = 0.11

# Conf: Compare intervals
pval[5] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 1], HR$resp1[HR$dist==2 & HR$resp2 == 1], paired=T)$p.value
pval[6] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 1], HR$resp1[HR$dist==3 & HR$resp2 == 1], paired=T)$p.value
pval[7] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 1], HR$resp1[HR$dist==4 & HR$resp2 == 1], paired=T)$p.value

pval[8] <- t.test(HR$resp1[HR$dist==2 & HR$resp2 == 1], HR$resp1[HR$dist==3 & HR$resp2 == 1], paired=T)$p.value # p = 0.005
pval[9] <- t.test(HR$resp1[HR$dist==2 & HR$resp2 == 1], HR$resp1[HR$dist==4 & HR$resp2 == 1], paired=T)$p.value # p = 0.004

pval[10] <- t.test(HR$resp1[HR$dist==3 & HR$resp2 == 1], HR$resp1[HR$dist==4 & HR$resp2 == 1], paired=T)$p.value

# Unconf: Compare intervals
pval[11] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 0], HR$resp1[HR$dist==2 & HR$resp2 == 0], paired=T)$p.value
pval[12] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 0], HR$resp1[HR$dist==3 & HR$resp2 == 0], paired=T)$p.value
pval[13] <- t.test(HR$resp1[HR$dist==1 & HR$resp2 == 0], HR$resp1[HR$dist==4 & HR$resp2 == 0], paired=T)$p.value

pval[14] <- t.test(HR$resp1[HR$dist==2 & HR$resp2 == 0], HR$resp1[HR$dist==3 & HR$resp2 == 0], paired=T)$p.value
pval[15] <- t.test(HR$resp1[HR$dist==2 & HR$resp2 == 0], HR$resp1[HR$dist==4 & HR$resp2 == 0], paired=T)$p.value

pval[16] <- t.test(HR$resp1[HR$dist==3 & HR$resp2 == 0], HR$resp1[HR$dist==4 & HR$resp2 == 0], paired=T)$p.value

pval_adj <- p.adjust(pval,'fdr')

# ANOVA across bins -------------------------------------------------------

# Load ezANOVA package
library("ez")

# Make bins to factors
data$dist=as.factor(data$dist)
data$ID=as.factor(data$ID)
data$resp1=as.factor(data$resp1)
data$resp2=as.factor(data$resp2)

# Call to ANOVA wrapper
output_anova = ezANOVA(data = HR,       # dataframe containing all relevant variables, see below
                       dv = resp1,        # dependent variable: detection (yes/no)
                       wid = ID,          # array with the participant ID within the dataframe
                       within= .(dist,resp2),   # specify the names of within-subject factors
                       type = 3,          # sum-of-squares-type, should be '3' in your case, but check function help to be sure
                       detailed = TRUE,   # some output options
                       return_aov = TRUE
)

print(output_anova)

# Perform test with corrected degrees of freedom
fit_all <- aov_ez("ID","resp1", HR, within=c("dist", "resp2"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values

fit_all$anova_table$`num Df`
fit_all$anova_table$`den Df`
fit_all$anova_table$`Pr(>F)`
fit_all$anova_table$F

# Main effect bin (dist) significant: F(3,120) = 6.98, p = 0.0003
# But Mauchly's Test for Sphericity significant
# Report Greenhouse-Geisser correction?
# df * GGe (Greenhouse-Geisser epsilon)
output_anova$ANOVA$DFn[4] * output_anova$`Sphericity Corrections`$GGe[2]
output_anova$ANOVA$DFd[4] * output_anova$`Sphericity Corrections`$GGe[2]
output_anova$ANOVA$F[2]
output_anova$`Sphericity Corrections`$`p[GG]`
# F(2.74,96.01) = 7.01, p = 0.0004


# Confidence across the cardiac cycle -------------------------------------

HR=aggregate(resp2 ~ ID*dist, data_select, FUN=mean)
HR=aggregate(resp2 ~ ID*dist*resp1, data_select, FUN=mean)

HR$ID <- as.factor(HR$ID)
HR$dist <- as.factor(HR$dist)
HR$resp1 <- as.factor(HR$resp1)


# Plot confidence rate across cardiac cycle -------------------------------

library(ggplot2)

# Boxplot of the hit rates for each time bin 
ggplot(HR, aes(x=dist, y=resp2, fill=resp1)) + # col=resp2
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Confidence rate x detection x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("Confidence rate") +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position=c(0.08,0.08)
  ) +
  scale_x_discrete(labels=c("1" = "0-200", "2" = "200-400", "3" = "400-600", "4" = "600-900")) +
  scale_fill_discrete(labels=c("miss", "hit")) +
  annotate(geom = "text", x = 1, y = 1.06, label = "p = 0.02", color = "black") +  # T1: conf vs unconf
  annotate("segment", x = 0.8, xend = 1.2, y = 1.03, yend = 1.03, colour = "black") +
  annotate(geom = "text", x = 2, y = 1.06, label = "p = 0.02", color = "black") +  # T2: conf vs unconf
  annotate("segment", x = 1.8, xend = 2.2, y = 1.03, yend = 1.03, colour = "black") +
  annotate(geom = "text", x = 3.2, y = 0.09, label = "p = 0.02", color = "black") +  # yes: T2 vs T4
  annotate("segment", x = 2.2, xend = 4.2, y = 0.06, yend = 0.06, colour = "black")

#  annotate(geom = "text", x = 3.7, y = -0.01, label = "p = 0.02", color = "black") +  # yes: T3 vs T4
#  annotate("segment", x = 3.2, xend = 4.2, y = -0.04, yend = -0.04, colour = "black") +

#  annotate(geom = "text", x = 2.7, y = 0.12, label = "p = 0.02", color = "black") + # yes: T1 vs T4
#  annotate("segment", x = 1.2, xend = 4.2, y = 0.09, yend = 0.09, colour = "black") +

  #annotate(geom = "text", x = 3, y = 1.06, label = "p = 0.04", color = "black") +  # T3: conf vs unconf
  #annotate("segment", x = 2.8, xend = 3.2, y = 1.03, yend = 1.03, colour = "black")

#+  stat_summary(aes(y = resp2,group=2), fun.y=mean, colour="red", geom="line",group=2)


# t-Tests confidence between bins -----------------------------------------

HR=aggregate(resp2 ~ ID*dist, data, FUN=mean)

t.test(HR$resp2[HR$dist==1], HR$resp2[HR$dist==2], paired=T)
t.test(HR$resp2[HR$dist==1], HR$resp2[HR$dist==3], paired=T)
t.test(HR$resp2[HR$dist==1], HR$resp2[HR$dist==4], paired=T) # p = 0.035

t.test(HR$resp2[HR$dist==2], HR$resp2[HR$dist==3], paired=T)
t.test(HR$resp2[HR$dist==2], HR$resp2[HR$dist==4], paired=T) # p = 0.012

t.test(HR$resp2[HR$dist==3], HR$resp2[HR$dist==4], paired=T) # p = 0.016

# Splitting by yes/no

HR=aggregate(resp2 ~ ID*dist*resp1, data, FUN=mean)


# Post-hoc t-tests: Confidence rate x detection x time since R peak -------

pval <- numeric(16)

# Yes vs. no confidence rate within interval
pval[1] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 0], HR$resp2[HR$dist==1 & HR$resp1 == 1], paired=T)$p.value # p = 0.003
pval[2] <- t.test(HR$resp2[HR$dist==2 & HR$resp1 == 0], HR$resp2[HR$dist==2 & HR$resp1 == 1], paired=T)$p.value # p = 0.002
pval[3] <- t.test(HR$resp2[HR$dist==3 & HR$resp1 == 0], HR$resp2[HR$dist==3 & HR$resp1 == 1], paired=T)$p.value # p = 0.053
pval[4] <- t.test(HR$resp2[HR$dist==3 & HR$resp1 == 0], HR$resp2[HR$dist==3 & HR$resp1 == 1], paired=T)$p.value # p = 0.22
#t.test(HR$resp2[HR$dist==4 & HR$resp1 == 0 & HR$ID != 6], HR$resp2[HR$dist==4 & HR$resp1 == 1 & HR$ID != 6], paired=T)$p.value # p = 0.22 / 0.23 N=40 ID06 no such case

# Yes: confidence rate between intervals
pval[5] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 1], HR$resp2[HR$dist==2 & HR$resp1 == 1], paired=T)$p.value
pval[6] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 1], HR$resp2[HR$dist==3 & HR$resp1 == 1], paired=T)$p.value
pval[7] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 1], HR$resp2[HR$dist==4 & HR$resp1 == 1], paired=T)$p.value # p = 0.02

pval[8] <- t.test(HR$resp2[HR$dist==2 & HR$resp1 == 1], HR$resp2[HR$dist==3 & HR$resp1 == 1], paired=T)$p.value
pval[9] <- t.test(HR$resp2[HR$dist==2 & HR$resp1 == 1], HR$resp2[HR$dist==4 & HR$resp1 == 1], paired=T)$p.value # p = 0.004

pval[10] <- t.test(HR$resp2[HR$dist==3 & HR$resp1 == 1], HR$resp2[HR$dist==4 & HR$resp1 == 1], paired=T)$p.value # p = 0.02

# No: confidence rate between intervals
pval[11] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 0], HR$resp2[HR$dist==2 & HR$resp1 == 0], paired=T)$p.value
pval[12] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 0], HR$resp2[HR$dist==3 & HR$resp1 == 0], paired=T)$p.value # p = 0.053
pval[13] <- t.test(HR$resp2[HR$dist==1 & HR$resp1 == 0], HR$resp2[HR$dist==4 & HR$resp1 == 0], paired=T)$p.value

pval[14] <- t.test(HR$resp2[HR$dist==2 & HR$resp1 == 0], HR$resp2[HR$dist==3 & HR$resp1 == 0], paired=T)$p.value
pval[15] <- t.test(HR$resp2[HR$dist==2 & HR$resp1 == 0], HR$resp2[HR$dist==4 & HR$resp1 == 0], paired=T)$p.value

pval[16] <- t.test(HR$resp2[HR$dist==3 & HR$resp1 == 0], HR$resp2[HR$dist==4 & HR$resp1 == 0], paired=T)$p.value

pval_adj <- p.adjust(pval,'fdr')

which(pval_adj<=0.05)

# ANOVA Confidence x bin --------------------------------------------------

library("ez")

HR=aggregate(resp2 ~ ID*dist*resp1, data, FUN=mean)

HR$ID <- as.factor(HR$ID)
HR$dist <- as.factor(HR$dist)
HR$resp1 <- as.factor(HR$resp1)

# Call to ANOVA wrapper
output_anova = ezANOVA(data = HR,       # dataframe containing all relevant variables, see below
                       dv = resp2,        # dependent variable: detection (yes/no)
                       wid = ID,          # array with the participant ID within the dataframe
                       within= .(dist,resp1),   # specify the names of within-subject factors
                       type = 3,          # sum-of-squares-type, should be '3' in your case, but check function help to be sure
                       detailed = TRUE,   # some output options
                       return_aov = TRUE
)

print(output_anova)

fit_all <- aov_ez("ID","resp2", HR, within=c("dist", "resp1"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values


# ANOVA -------------------------------------------------------------------

HR_wo6 <- subset(HR, ID != 6)
HR_wo6$ID <- as.factor(HR_wo6$ID)

HR <- HR_wo6

library("afex")

fit_all <- aov_ez("ID","resp1", HR, within="dist")
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values

# Plot hit rates across bins -----------------------------------------------

library(ggplot2)
library("rstatix")
library("dplyr")
library("ggpubr")

library("tidyverse")

# Create labels based on intervals
label_list <- paste(as.character(head(intervals,-1)*1000),as.character(intervals[-1]*1000),sep = "\n-")
names(label_list) <- 1:(length(intervals)-1)


# Test PLV
stat_test = HR %>% 
  t_test(resp1 ~ dist, paired=TRUE, p.adjust.method = "fdr", detailed=TRUE) %>%
  add_xy_position(step.increase=.03)
  add_xy_position(dodge = 0.8)
  

# Boxplot of the hit rates for each time bin 
ggplot(HR, aes(x=dist, y=resp1)) +
  geom_boxplot(aes(x=dist, y=resp1, fill=rgb(190,186,218, maxColorValue=255) )) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Detection rate x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("Detection rate") +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.position = "none"
    #legend.title = element_blank(),
    #legend.text = element_text(size = 10),
    #legend.position=c(0.12,0.10)
  ) +
  stat_pvalue_manual(stat_test, label = "p.adj", tip.length = 0.01, hide.ns=T, step.increase = 0.1) +
  scale_x_discrete(labels=label_list)


# Mean curve --------------------------------------------------------------

HR_mean <- aggregate(resp1 ~ dist, HR, FUN=mean)
HR_mean$dist <- as.numeric(HR_mean$dist)    

# Create labels based on intervals
label_list <- paste(as.character(head(intervals,-1)*1000),as.character(intervals[-1]*1000),sep = "\n-")
names(label_list) <- 1:(length(intervals)-1)

ggplot(HR_mean, aes(dist,resp1)) + geom_point() + geom_smooth() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(1,length(intervals)-1,1), labels = label_list) +
  ggtitle("Detection rate x time since R peak") +
  xlab("Time since R-peak in ms") + ylab("Detection rate") +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
  ) +
  geom_hline(yintercept = .5, size = 0.3)
