# respirationCA - Plot Response Times by Confidence -----------------------------------------

# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      June 2, 2020

# DATA - resp1 --------------------------------------------------------------------

# Mean response times
#null_resp1_t_no_conf = dID1$null_resp1_t_no_conf,
#null_resp1_t_no_unconf = dID1$null_resp1_t_no_unconf,
resp1_t <- data.frame(near_resp1_t_no_conf = dID$near_resp1_t_no_conf,
                      near_resp1_t_no_unconf = dID$near_resp1_t_no_unconf,
                      near_resp1_t_yes_conf = dID$near_resp1_t_yes_conf,
                      near_resp1_t_yes_unconf = dID$near_resp1_t_yes_unconf)

# Number of trials
#null_n_no_conf = dID1$null_n_no_conf*dID1$block_num,
#null_n_no_unconf = dID1$null_n_no_unconf*dID1$block_num,
trial_n <- data.frame(near_n_no_conf = dID$near_n_no_conf,
                      near_n_no_unconf = dID$near_n_no_unconf,
                      near_n_yes_conf = dID$near_n_yes_conf,
                      near_n_yes_unconf = dID$near_n_yes_unconf)

# Title
main_title <- substitute(paste('Yes/no response time (N = ', x, ')', sep = ''), list(x=nrow(resp1_t)))
main_title <- 'Yes/no response time'

# Category labels
#cond_lab <- c('Miss confident', 'Miss unconfident', 'Hit confident', 'Hit unconfident')
cond_lab <- c('Miss C+', 'Miss C-', 'Hit C+', 'Hit C-')

# DATA - resp2 --------------------------------------------------------------------

# Mean response times
resp1_t <- data.frame(near_resp2_t_no_conf = dID$near_resp2_t_no_conf,
                      near_resp2_t_no_unconf = dID$near_resp2_t_no_unconf,
                      near_resp2_t_yes_conf = dID$near_resp2_t_yes_conf,
                      near_resp2_t_yes_unconf = dID$near_resp2_t_yes_unconf)

# Number of trials
trial_n <- data.frame(near_n_no_conf = dID$near_n_no_conf,
                      near_n_no_unconf = dID$near_n_no_unconf,
                      near_n_yes_conf = dID$near_n_yes_conf,
                      near_n_yes_unconf = dID$near_n_yes_unconf)

# Title
main_title <- substitute(paste('Confidence response time (N = ', x, ')', sep = ''), list(x=nrow(resp1_t)))
main_title <- 'Confidence response time'

# Category labels
#cond_lab <- c('Miss confident', 'Miss unconfident', 'Hit confident', 'Hit unconfident')
cond_lab <- c('Miss C+', 'Miss C-', 'Hit C+', 'Hit C-')


# SETTINGS ----------------------------------------------------------------

# Y-axis label
y_axis_label <- 'Mean response time (s)'

# Category labels (x-axis)
x_labels <- c(paste(cond_lab[1],'\n(n~', round(mean(trial_n$near_n_no_conf)), ')', sep = ''),
              paste(cond_lab[2],'\n(n~', round(mean(trial_n$near_n_no_unconf)), ')', sep = ''),
              paste(cond_lab[3],'\n(n~', round(mean(trial_n$near_n_yes_conf)), ')', sep = ''),
              paste(cond_lab[4],'\n(n~', round(mean(trial_n$near_n_yes_unconf)), ')', sep = ''))
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
        col=colmap[c(3,3,4,4)],
        range = 1.5,
        ylim = c(0,1.1),
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
     at = 1:4,
     labels = x_labels)#,
#mgp = c(3, 3, 0))


# PLOT T-TEST - resp1 -------------------------------------------------------------

# Miss confident vs. miss unconfident
t.test(dID$near_resp1_t_no_conf-dID$near_resp1_t_no_unconf)
# t = -15.462, df = 40, p-value < 2.2e-16

lines(c(1,2),c(0.30,0.30), col = 1, lwd = 1)

text(1.5,0.29,'p < 0.001', pos = 3, col = 1, cex = 1)

# Hit confident vs. hit unconfident
t.test(dID$near_resp1_t_yes_conf-dID$near_resp1_t_yes_unconf)
# t = -11.799, df = 40, p-value = 1.329e-14

lines(c(3,4),c(0.22,0.22), col = 1, lwd = 1)

text(3.5,0.21,'p < 0.001', pos = 3, col = 1, cex = 1)


# PLOT T-TEST - resp2 -------------------------------------------------------------

# Miss confident vs. miss unconfident
t.test(dID$near_resp2_t_no_conf-dID$near_resp2_t_no_unconf)
# t = -5.8901, df = 40, p-value = 6.747e-07

lines(c(1,2),c(0.05,0.05), col = 1, lwd = 1)

#text(1.5,0.06,'p = 0.0001', pos = 3, col = 1, cex = 1)
text(1.5,0.04,'p < 0.001', pos = 3, col = 1, cex = 1)

# Hit confident vs. hit unconfident
t.test(dID$near_resp2_t_yes_conf-dID$near_resp2_t_yes_unconf)
# t = -1.3044, df = 40, p-value = 0.1995

