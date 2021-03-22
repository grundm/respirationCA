# respirationCA - Plot Response Times -----------------------------------------

# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      June 2, 2020

# DATA - resp1 --------------------------------------------------------------------

# Mean response times
resp1_t <- data.frame(null_resp1_t_no = dID$null_resp1_t_no,
                      near_resp1_t_no = dID$near_resp1_t_no,
                      near_resp1_t_yes = dID$near_resp1_t_yes)

# Number of trials
trial_n <- data.frame(null_n_no = dID$null_n_no,
                      near_n_no = dID$near_n_no,
                      near_n_yes = dID$near_n_yes)

# Title
main_title <- substitute(paste('Yes/no response time (N = ', x, ')', sep = ''), list(x=nrow(resp1_t)))
main_title <- 'Yes/no response time'

# Category labels
cond_lab <- c('CR', 'Miss', 'Hit')

# DATA - resp2 --------------------------------------------------------------------

# Mean response times
resp1_t <- data.frame(null_resp1_t_no = dID$null_resp2_t_no,
                      near_resp1_t_no = dID$near_resp2_t_no,
                      near_resp1_t_yes = dID$near_resp2_t_yes)

# Number of trials
trial_n <- data.frame(null_n_no = dID$null_n_no,
                      near_n_no = dID$near_n_no,
                      near_n_yes = dID$near_n_yes)

# Title
main_title <- substitute(paste('Confidence response time (N = ', x, ')', sep = ''), list(x=nrow(resp1_t)))
main_title <- 'Confidence response time'

# Category labels
cond_lab <- c('CR', 'Miss', 'Hit')


# SETTINGS ----------------------------------------------------------------

# Y-axis label
y_axis_label <- 'Mean response time (s)'

# Category labels (x-axis)
x_labels <- c(paste(cond_lab[1],'\n(n~', round(mean(trial_n$null_n_no)), ')', sep = ''),
              paste(cond_lab[2],'\n(n~', round(mean(trial_n$near_n_no)), ')', sep = ''),
              paste(cond_lab[3],'\n(n~', round(mean(trial_n$near_n_yes)), ')', sep = ''))
x_labels <- cond_lab

# Create colormap
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=4
colmap <- c(rgb(141,211,199, maxColorValue=255),
            rgb(255,255,179, maxColorValue=255),
            rgb(190,186,218, maxColorValue=255),
            rgb(251,128,114, maxColorValue=255))


# PLOT DETECTION ----------------------------------------------------------

# Plot settings
par(cex.main = 2.4,
    cex.lab = 2.0,
    cex.axis = 1.8,
    #mar = c(4.5, 5.5, 4, 1),
    mar = c(2.5, 5.5, 4, 1),
    yaxs = 'i',
    xaxs = 'i',
    yaxt = 'n',
    xaxt = 'n',
    bty = 'o')

# Draw boxplot
boxplot(resp1_t,
        col=colmap[c(1,3,4)],
        range = 1.5,
        ylim = c(0,1),
        whisklty = 1,
        outpch = 21,
        outbg = 'black',
        outcex = .9)

# Set title
title(main = main_title)

# Label y-axis
title(ylab = y_axis_label,
      line = 3.6)

# Activitate axes
par(xaxt = 's')
par(yaxt = 's')

axis(side = 2,
     at = seq(0, 1, .2),
     las = 1)

axis(side = 1,
     at = 1:3,
     labels = x_labels)#,
     #mgp = c(3, 3, 0))


# PLOT T-TEST - resp1 -------------------------------------------------------------

# CR vs. miss
t.test(dID$null_resp1_t_no-dID$near_resp1_t_no)
# t = -6.8507, df = 40, p-value = 3.028e-08

lines(c(1,2),c(0.30,0.30), col = 1, lwd = 1)

text(1.5,0.29,'p < 0.001', pos = 3, col = 1, cex = 1)

# CR vs. hit
t.test(dID$null_resp1_t_no-dID$near_resp1_t_yes)
# t = -2.3934, df = 40, p-value = 0.02147

lines(c(1,3),c(0.22,0.22), col = 1, lwd = 1)

text(2,0.21,'p = 0.02', pos = 3, col = 1, cex = 1)

# Miss vs. hit
t.test(dID$near_resp1_t_no-dID$near_resp1_t_yes)
# t = 0.7972, df = 40, p-value = 0.43


t.test(dID$near_resp1_t_conf-dID$near_resp1_t_unconf)

data_stims = readRDS(paste(data_dir, "resp_data_stims_martin_20210208.Rds", sep="/"))

data_stims_valid <- subset(data_stims, resp_filter == 1 & HR_larger_FAR_filter == 1)

resp_times <- aggregate(resp1_t ~ trial_type*ID, subset(data_stims_valid, trial_type != "FA"), FUN = mean)

stat_test = resp_times %>% 
  t_test(resp1_t ~ trial_type, paired=TRUE, p.adjust.method = "fdr")

aggregate(resp1_t ~ trial_type, resp_times, FUN = mean)

# PLOT T-TEST - resp2 -------------------------------------------------------------

# CR vs. miss
t.test(dID$null_resp2_t_no-dID$near_resp2_t_no)
# t = -3.8503, df = 40, p-value = 0.000416

lines(c(1,2),c(0.07,0.07), col = 1, lwd = 1)

#text(1.5,0.06,'p = 0.0004', pos = 3, col = 1, cex = 1)
text(1.5,0.06,'p < 0.001', pos = 3, col = 1, cex = 1)

# CR vs. hit
t.test(dID$null_resp2_t_no-dID$near_resp2_t_yes)
# t = -0.17539, df = 40, p-value = 0.8617


# Miss vs. hit
t.test(dID$near_resp2_t_no-dID$near_resp2_t_yes)
# t = 1.5666, df = 40, p-value = 0.1251

