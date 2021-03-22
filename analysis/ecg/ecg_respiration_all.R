require(R.matlab)
rm(list = ls())
setwd('D:/respirationCA/ecg_analysis/kubios')
bigdata=c()
alldata=read.delim("D:/respirationCA/data/respirationCA_trials.txt")
subids=c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '1001')

#subject 16 remove heart interval with 1.6 s length, premature R
#subject 17 has immature heartbeats
#subject 23-block 1_ 1 immature heartbeat


#subject 31 block 2 and block 4--- no ecg signal
#subject 37 block 4, kubios did not recognize R peaks (negative ecg signal)
#so subjects 31 and 37 don't have kubios files 

#subject 3 block 2, first 4 trials don't have ecg recording
sub3del=which(alldata$ID==3 & alldata$block==2 & alldata$trial==1)
sub3del=sub3del:(sub3del+150-146-1)

#subject 21 block 1, some first trials don't have ecg recording
sub22del=which(alldata$ID==22 & alldata$block==1 & alldata$trial==1)
sub22del=sub22del:(sub22del+150-128-1)

alldata=alldata[-c(sub3del, sub22del),]

alldata$block[alldata$ID==21 & alldata$block==2]=1 ### subject 21 block 1 and block 2 are recorded together!, therefore kubios file labeled as Kubios_ID21_1 contains both block 1 and 2.

for (s in 1:length(subids)){
  subdata=c()
  for (b in 1:4){
    kubiosname=paste('Kubios_ID', subids[s],'_',b, '_hrv.mat', sep='')
    if (! file.exists(kubiosname)) next # if no kubios file is present, go to the next
    data=c()
    if (s==41){
      data=subset(alldata, ID==subids[s] & block==b)
    } else {
      data=subset(alldata, ID==s & block==b)
    }
    stims=readMat(paste('ecg',subids[s],'_', b, '.mat', sep=''))
    stims=stims$ecg.event
    data$stim_onset=t(stims)
    
    R_peaks= readMat(kubiosname)
    R_peaks= as.vector(R_peaks$Res[[4]][[2]][[2]]) 
    #initialize variables of the data frame
    data$stim_degree=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$diff2peak=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$RR_int_minus2=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$RR_int_minus1=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$RR_int_zero_stimulus=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$RR_int_plus1=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$RR_int_plus2=matrix(data=NA,nrow=length(data$trial),ncol=1)
    data$RR_int_plus3=matrix(data=NA,nrow=length(data$trial),ncol=1)
    for (ind in 1:length(data$trial)) { 
      # the time of the peak just before stimulus onset
      pos=max(which(R_peaks < data$stim_onset[ind]))
      # the difference between the onset and previous R peak
      diff2peak=data$stim_onset[ind] - R_peaks[pos]
      data$diff2peak[ind]=diff2peak
      # relative position of onset in cardiac cycle
      stim_degree=360 * diff2peak/(R_peaks[pos+1] - R_peaks[pos])
      data$stim_degree[ind]=stim_degree
      
      
      #heart rate calculations
      if (pos>2){data$RR_int_minus2[ind]=R_peaks[pos-1] - R_peaks[pos-2]}
      if (pos>1){data$RR_int_minus1[ind]=R_peaks[pos] - R_peaks[pos-1] }
      data$RR_int_zero_stimulus[ind]=R_peaks[pos+1] - R_peaks[pos] 
      data$RR_int_plus1[ind]=R_peaks[pos+2] - R_peaks[pos+1] 
      data$RR_int_plus2[ind]=R_peaks[pos+3] - R_peaks[pos+2] 
      data$RR_int_plus3[ind]=R_peaks[pos+4] - R_peaks[pos+3] 
    }
    subdata=rbind(subdata,data)
  }
  bigdata=rbind(bigdata,subdata)
}
write.csv2(bigdata,file="ecgdata_all.csv")  ###### I send you this file #######


