# respirationCA - Plot Detection Rates ------------------------------------------

# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      June 2, 2020

# DATA - ALL --------------------------------------------------------------------

# Mean yes response rate
yes_resp <- data.frame(null_yes = dID$null_resp1,
                       near_yes = dID$near_resp1)

# Number of trials
trial_n <- data.frame(null_n = dID$null_n,
                      near_n = dID$near_n)

# Title
main_title <- substitute(paste('Detection (N = ', x, ')', sep = ''), list(x=nrow(yes_resp)))

# Category labels
cond_lab <- c('Catch', 'Near')


# SETTINGS ----------------------------------------------------------------

# Y-axis label
y_axis_label <- 'P(yes)'

# Category labels (x-axis)
x_labels <- c(paste(cond_lab[1],'\n(n~', round(mean(trial_n$null_n)), ')', sep = ''),
              paste(cond_lab[2],'\n(n~', round(mean(trial_n$near_n)), ')', sep = ''))

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
boxplot(yes_resp,
        col=colmap[c(1,4)],
        range = 1.5,
        ylim = c(0,1),
        whisklty = 1,
        outpch = 21,
        outbg = 'black',
        outcex = .9)

# # Draw line for 50%
# abline(h = 0.5,
#        lty = 2,
#        col = 'black')
# 
# # Label line
# text(.5,.53,'50%', pos = 4, col = 'black', cex = 1.7)

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
     at = 1:2,
     labels = x_labels,
     mgp = c(3, 3, 0))