# respirationCA - Plot Confidence Rates -----------------------------------------

# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      June 2, 2020

# DATA - ALL --------------------------------------------------------------------

# Mean confidence rate
resp2 <- data.frame(null_resp2_no = dID$null_resp2_no,
                    near_resp2_no = dID$near_resp2_no,
                    near_resp2_yes = dID$near_resp2_yes)

# Number of trials
trial_n <- data.frame(null_n_no = dID$null_n_no,
                      near_n_no = dID$near_n_no,
                      near_n_yes = dID$near_n_yes)

# Title
main_title <- substitute(paste('Confidence (N = ', x, ')', sep = ''), list(x=nrow(resp2)))

# Category labels
cond_lab <- c('CR', 'Miss', 'Hit')


# SETTINGS ----------------------------------------------------------------

# Y-axis label
y_axis_label <- 'P(confident)' #'Mean rate of "confident"'

# Category labels (x-axis)
x_labels <- c(paste(cond_lab[1],'\n(n~', round(mean(trial_n$null_n_no)), ')', sep = ''),
              paste(cond_lab[2],'\n(n~', round(mean(trial_n$near_n_no)), ')', sep = ''),
              paste(cond_lab[3],'\n(n~', round(mean(trial_n$near_n_yes)), ')', sep = ''))

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
    mar = c(4.5, 5.5, 4, 1),
    yaxs = 'i',
    xaxs = 'i',
    yaxt = 'n',
    xaxt = 'n',
    bty = 'o')

# Draw boxplot
boxplot(resp2,
        col=colmap[c(1,3,4)],
        range = 1.5,
        ylim = c(0,1),
        whisklty = 1,
        outpch = 21,
        #outcol = 'darkgrey',
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
     labels = x_labels,
     mgp = c(3, 3, 0))


# PLOT T-TEST -------------------------------------------------------------

# CR vs. miss
t.test(dID$null_resp2_no-dID$near_resp2_no)
# t = 6.8478, df = 40, p-value = 3.057e-08

lines(c(1,2),c(0.15,0.15), col = 1, lwd = 1)

text(1.5,0.14,'p < .001', pos = 3, col = 1, cex = 1)
#text(1.5,0.13,'p = 0.00000003', pos = 3, col = 1, cex = 1)

# CR vs. hit
t.test(dID$null_resp2_no-dID$near_resp2_yes)
# t = 7.9009, df = 40, p-value = 1.078e-09

lines(c(1,3),c(0.04,0.04), col = 1, lwd = 1)

text(2,0.03,'p < .001', pos = 3, col = 1, cex = 1)
#text(2,0.06,'p = 0.000000001', pos = 3, col = 1, cex = 1)

# Miss vs. hit
t.test(dID$near_resp2_no-dID$near_resp2_yes)
# t = 2.4378, df = 40, p-value = 0.01932

lines(c(2,3),c(0.11,0.11), col = 1, lwd = 1)

text(2.5,0.10,'p = 0.019', pos = 3, col = 1, cex = 1)
