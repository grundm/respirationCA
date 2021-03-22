# respirationCA - Stimulus onset relative to R peak -----------------------

# Author:         Esra Al, Martin Grund
# Last update:    November 6, 2019

# Loads the R peaks (preprocessed ECG data with Kubios) and puts it in relation
# to the stimulus onsets (trigger timestamps). Calculates the time difference
# of a stimulus onset to the previous R peak and its degree on the cardiac
# cycle.

# Settings ----------------------------------------------------------------

rm(list = ls())

# Behavioral data directory
data_dir <- '~/ownCloud/promotion/experiment/respirationCA/data'

# Kubios data directory
kubios_path <- 'D:/respirationCA/ecg_analysis/kubios'

# Behavioral data
behav_data <- paste(data_dir, '/respirationCA_trials.txt', sep = '')

# Output
output_file <- paste(data_dir, '/ecgdata.csv', sep = '')

# Define participant IDs
subids=c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
         '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
         '21', '22', '23', '24', '25', '26', '27', '28', '29', '30',
         '31', '32', '33', '34', '35', '36', '37', '38', '39', '40',
         '1001')

require(R.matlab)


# Load data ---------------------------------------------------------------

alldata <- read.delim(behav_data)


# Quality assessment ------------------------------------------------------

# subject 16 remove heart interval with 1.6 s length, premature R
# subject 17 has immature heartbeats
# subject 23-block 1_ 1 immature heartbeat

# subject 31 block 2 and block 4--- no ecg signal
# subject 37 block 4, kubios did not recognize R peaks (negative ecg signal)

# subject 3 block 2, first 4 trials don't have ecg recording
sub3del=which(alldata$ID==3 & alldata$block==2 & alldata$trial==1)
sub3del=sub3del:(sub3del+150-146-1)

# subject 22 block 1, some first trials don't have ecg recording
sub22del=which(alldata$ID==22 & alldata$block==1 & alldata$trial==1)
sub22del=sub22del:(sub22del+150-128-1)

# Removes trials with missing triggers for ID03 and ID22
alldata=alldata[-c(sub3del, sub22del),]

# subject 21 block 1 and block 2 are recorded together!!!!!!
# -> Set block == 1?
alldata$block[alldata$ID==21 & alldata$block==2]=1

# It's a bit confusing that the behavioral data is deleted or edited. Necessary?


# Calculate time to previous R peak and degree on cardiac cycle ------------

# Prepare output variable
bigdata=c()

# Loop participants
for (s in 1:length(subids)){
  
  # Prepare variable for participants output data?
  subdata=c()
  
  # Loop blocks
  for (b in 1:4){
    
    # Kubios output file (preprocessed ECG data = R peaks)
    kubiosname = paste(kubios_path, '/Kubios_ID', subids[s],'_',b, '_hrv.mat', sep='')
    
    # Check if file exists
    if (! file.exists(kubiosname)) next
    data=c()
    # Select behavioral data from current participants s and current block b
    data=subset(alldata, ID==subids[s] & block==b, select = c('ID','block','trial','stim_type','resp1', 'resp2'))
    
    # Read ECG data
    stims = readMat(paste(kubios_path, '/ecg',subids[s],'_', b, '.mat', sep=''))
    # Get the trigger onsets
    stims=stims$ecg.event
    data$stim_onset=t(stims)
    
    # Read Kubios file?
    R_peaks= readMat(kubiosname)
    # Get R peaks
    R_peaks= as.vector(R_peaks$Res[[4]][[2]][[2]])
    
    # Prepare variables for degree and difference to previous R peak
    data$stim_degree=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$diff2peak=matrix(data=NA,nrow=length(data$trial),ncol=1)
    
    # Loop trials
    for (ind in 1:length(data$trial)) { 
      
      # Get time of the R peak before the stimulus onset
      pos=max(which(R_peaks < data$stim_onset[ind]))
      
      # Calculate difference between the stimulus onset and the previous R peak
      diff2peak=data$stim_onset[ind] - R_peaks[pos]
      data$diff2peak[ind]=diff2peak
      
      # Calculate relative position of the stimulus onset on the cardiac cycle
      stim_degree=360 * diff2peak/(R_peaks[pos+1] - R_peaks[pos])
      data$stim_degree[ind]=stim_degree
    }
    
    # Append current block data to current participant's data
    subdata=rbind(subdata,data)
  }
  # Append current participant's data to all participants
  bigdata=rbind(bigdata,subdata)
}

# Write the data to the file
write.csv2(bigdata, file = output_file)