# respirationCA - t-wave detection ----------------------------------------


# Author:         Martin Grund, Esra Al
# Last update:    June 28, 2021

# Adapted script: https://github.com/Esra-Al/systole_detection/blob/main/twave_algorithm.R

detect_twave <- function(s,b,trial,R_peak,RR_int_s,ECG,ecg_srate,print_jpeg,file_suffix) {

# Input variables:    
# s = participant index
# b = block number
# trial - trial number
# R_peak in s
# R_int in s
# ECG signal as two row vector (first row: time in s; second row: signal)
# ecg_srate in samples/s (Hz)
# print_jpeg = TRUE/FALSE
# file_suffix - string to indicate RR-interval
  
# T-wave detection settings --------------------------------------------------------

# Cut-off end of RR-interval to exclude next R-peak in s
# as factor of sampling frequency (e.g., 5000 Hz * 0.06 s = 300 samples)
cut_off_factor <- 0.06 # = 60 ms

# R interval is extended for visualization purposes to include both R-peaks by 120 ms
# as factor of sampling frequency (e.g., 5000 Hz * 0.12 s = 600 samples)
twave1_padding <- c(0.12, 0.12)

# When t-wave should occur latest in s
twave_onset_max <- 0.14
if (s == 6) { twave_onset_max <- 0.07 }

# Cut-off RR-interval by 140 ms and devide remaining vector
twave_peak_cutoff <- 0.14
twave_peak_interval_factor <- 3

# T-wave end search window from t-max in s
twave_end_search <- 0.14
twave_inflection_search <- 0.14

# Time window from xm to xr (xm = inflection point point after T-peak in twave_inflection_search window (minimum value in first derivative)
twave_end_window <- 0.08
if (s == 6) { twave_end_window <- 0.04 }

# Function: find t-wave end -----------------------------------------------

# First, calculation of the trapeziums areas of all the points located between xseq and yseq,
# then identification of the point which gives the maximum area and label it as the t-wave end.

trapez_area <- function(xm, ym, xseq, yseq, xr) {
  a <- numeric()
  for (i in seq_along(xseq)){
    a[i] <- 0.5 * (ym - yseq[i]) * ((2*xr) - xseq[i] - xm)
  }
  x_tend <- which.max(a)+xm-1
  return(x_tend)
}


# Determine t-wave and hence systole and diastole duration ----------------

# Calculate position of R peak that occurred just before the stimulus onset in data points
ecgpos1 = R_peak*ecg_srate #fs=1000/5000

# Calculate the duration of RR interval containing the stimulus in data points
RRint = RR_int_s*ecg_srate #in data points

# Calculate a time interval where the end of t wave would be included (in data points).
# Cut off fs*(150/2500) data points so that t wave max calculation would not detect the next R peak.    
ecgpos2 = ecgpos1 + RRint - (ecg_srate*cut_off_factor)

# create a vector containing all possible data points in the RR interval for visualization purposes.
# 120 ms or fs*(300/2500) data points are chosen arbitrarily to include the Rpeaks in the plots.
twave1 = ECG[(ecgpos1-(ecg_srate*twave1_padding[1])):(ecgpos2+(ecg_srate*twave1_padding[2])), ]
# Visual check: plot(twave1,type='l')

# Create a vector (twave) where the maximum of t wave could be calculated. 
# After visually observing a block of each participant, it is decided that tmax should occur at least after 350/2500 ms after the R peak.
twave = ECG[(ecgpos1+(ecg_srate*twave_onset_max)):ecgpos2, ]
# Visual check: plot(twave,type='l')

# In order to calculate the position of tmax, delete (350/2500)*fs data points from the end of RRint vector 
# and search it until the 1/3 of the remaining vector.
# 1/3 ratio was given to compensate for the variation of trial to trial RR interval durations.

# Calculate a snip of t-wave which starts with the t-max
tmaxpos = which.max(twave[1:( (RRint-(ecg_srate*twave_peak_cutoff)) / twave_peak_interval_factor),2])
# Visual check: plot(twave[1:( (RRint-(ecg_srate*twave_peak_cutoff)) / twave_peak_interval_factor), ], type = 'l')
twave2 = twave[tmaxpos:dim(twave)[1],]
# Visual check: plot(twave2,type='l')


# Determine a point called xm located in the segment after the T peak, which has a minimum value in the first derivative.
# The algorithm searches for xm in a 120 ms time window starting from tmax.
# In case twave2 does not contain 0.12*fs data points, it searches only until the last point of twave2.
dp = twave_inflection_search*ecg_srate

if (dp>dim(twave2)[1]) {
  xm=which(diff(twave2[,2])==min(diff(twave2[,2])))
} else {
  xm=which(diff(twave2[1:dp,2])==min(diff(twave2[1:dp,2])))
}
xm=xm[1]
ym=twave2[xm,2]


# Determine a point xr which is supposed to happen after tend.
xr = xm + twave_end_window*ecg_srate

# Vector starting from xm and goes until xr
xseq=xm:xr
# Skip trial if t-wave was not properly detected
if (xr > dim(twave2)[1]) return(data.frame(twave_end=NA,syslen=NA,dias_start=NA,diaslen=NA))

yseq=twave2[xm:xr,2]

# Search t-wave end
tend = trapez_area(xm, ym, xseq, yseq, xr)
twave_end = twave2[tend,1]

# Next R-peak - 4 ms
nextR = (R_peak+RR_int_s-0.004)*ecg_srate

# Systole length (R-peak until t-wave end)
syslen = twave2[tend,1]-R_peak

# Diastole time sample vector and ecg trace
#diastole_dp = (nextR-(syslen*ecg_srate)):(nextR-1)
diastole_dp = ((twave2[tend,1]*ecg_srate)+1):(nextR-1)
diastole_ecg = ECG[diastole_dp, ]

dias_start = diastole_ecg[1,1]
diaslen = RR_int_s - syslen

#To check whether twave end detection worked well

if (print_jpeg) { jpeg(file = paste(kubios_dir, '/twave_detection/subject_', s,'_block_',b,'_trial_', trial, '_', file_suffix, '.jpg', sep='')) }
par(mfrow=c(1,2))
plot(twave1,col='black',xlab='time(ms)', ylab='electrical potential(mV)', type = 'l')
points(twave[tmaxpos,1],twave[tmaxpos,2],col='magenta',pch='+',cex=2) # t-wave peak
points(twave2[tend,1],twave2[tend,2],col='green',pch='+',cex=2) # t-wave end
points(diastole_ecg[1,1],diastole_ecg[1,2],col='red',pch='+',cex=2) # distole start
points(diastole_ecg[1,1],diastole_ecg[dim(diastole_ecg)[1],2],col='brown',pch='+',cex=2) # distole start
points(diastole_ecg[dim(diastole_ecg)[1],1],diastole_ecg[dim(diastole_ecg)[1],2],col='brown',pch='+',cex=2) # distole end


plot(twave2,col='black', xlab='time(ms)', ylab='electrical potential(mV)', type = 'l')
title(paste('subject ',s,' block ',b,' trial ', trial, ' RR ', file_suffix, sep=''),line=-2, outer=TRUE)
points(twave2[xm,1],twave2[xm,2],col='blue',pch='+',cex=2)
points(twave2[xr,1],twave2[xr,2],col='blue',pch='+',cex=2)
points(twave2[tend,1],twave2[tend,2],col='green',pch='+',cex=2) # t-wave end
points(twave[tmaxpos,1],twave[tmaxpos,2],col='magenta',pch='+',cex=2) # t-wave peak
#title(paste('subject ',s,'_block ',b,'_trial ', trial, 'RR_', file_suffix, sep=''),line=-2, outer=TRUE)

if (print_jpeg) { dev.off() }


return(data.frame(twave_end,syslen,dias_start,diaslen))

}