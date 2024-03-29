plot_test_R <- function (data,stim_type,resp1,resp2,title_str) {

  library(circular)
    
  # Get IDs
  ID_list = unique(data$ID)
  
  # Average data 
  mean_angle = matrix(data=NA,nrow=length(ID_list),ncol=1)
  var_angle = matrix(data=NA,nrow=length(ID_list),ncol=1)
  mean_length = matrix(data=NA,nrow=length(ID_list),ncol=1)
  r_pval = matrix(data=NA,nrow=length(ID_list),ncol=1)
  HR = matrix(data=NA,nrow=length(ID_list),ncol=1)
  
  for (s in 1:length(ID_list)) {
    angles_tmp <- circular(data$stim_degree[data$ID==ID_list[s] & data$stim_type %in% stim_type & data$resp1 %in% resp1 & data$resp2 %in% resp2], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0)
    #angles_tmp <- circular(data$stim_degree[data$ID==ID_list[s] & data$stim_type %in% stim_type], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0)
    mean_angle[s] = mean(angles_tmp, na.rm = TRUE)
    var_angle[s] = var(angles_tmp, na.rm = TRUE)
    #mean_angle[s] = mean(circular(data$stim_degree[data$ID==ID_list[s] & data$stim_type %in% stim_type & data$resp1 %in% resp1 & data$resp2 %in% resp2], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0), na.rm = TRUE)
    #var_angle[s] = var(circular(data$stim_degree[data$ID==ID_list[s] & data$stim_type %in% stim_type & data$resp1 %in% resp1 & data$resp2 %in% resp2], type="angles", units="degree", rotation="clock", modulo = "2pi", zero=0), na.rm = TRUE)
    mean_length[s] = rho.circular(angles_tmp, na.rm=T);
    r_test_tmp <- rayleigh.test(angles_tmp);
    r_pval[s] <- r_test_tmp$p.value;
    HR[s] <- mean(data$resp1[data$ID==ID_list[s] & data$stim_type==1])
  }
  
  print('Number of participants with significant Rayleigh test:')
  print(sum(r_pval < 0.05))
  print('FDR-corrected:')
  print(sum(p.adjust(r_pval,'fdr') < 0.05))
  
  
  # Plot mean angle distribution
  angle  =  circular(mean_angle, type="angles", units="degree", rotation="clock", modulo = "2pi", zero=pi/2)
  
  #par(mar = c(2, 2, 2, 2))
  par(mar = c(2, 2, 6, 2))
  # Inhale onsets:
  #par(mar = c(6, 2, 4, 2))
  plot(angle,
       stack=TRUE,
       bins = 720,
       col = "gray35",
       cex = 1.6)#,
       #shrink = 1.05)
  
  #title(title_str, line = 0)
  title(title_str, line = 4, cex.main = 2)
  # Inhale onsets:
  #title(title_str, line = 2, cex.main = 2)
  
  # Arrow 
  #arrows.circular(angle, col = "grey80", lwd = 2)
  arrows.circular(angle, y = mean_length, col = "grey80", lwd = 2, code = 0)
  arrows.circular(mean(angle,na.rm=T), y=rho.circular(angle, na.rm=T), lwd=5, col = "grey25")
  circ.dens = density(angle[!is.na(angle)], bw=20)
  #lines(circ.dens, col="grey25", lwd = 3, xpd=TRUE, shrink = 2) 
  lines(circ.dens, col="grey25", lwd = 4, xpd=TRUE) 
  
  # Print number of participants without a value
  print('Number of participants with NA:')
  print(sum(is.na(angle)))
  
  # Print number of participants with mean angle a-b
  print('Number of participants with mean angle in last quarter (270-0°):')
  #print(angle)
  #print(sum(-90 < angle & angle < 0, na.rm = TRUE))
  print(sum(angle > 270 & angle < 360, na.rm = TRUE))
  
  print('Number of participants with mean angle in 2nd quarter (90-180°):')
  #print(angle)
  print(sum(angle > 90  & angle < 180, na.rm = TRUE))
  
  print('Mean angle:')
  print(mean(angle, na.rm = TRUE))
  
  # Rayleigh test for uniformity
  r_test <- rayleigh.test(angle)
  print(r_test)
  
  text(0,-0.45,paste('R = ', round(r_test$statistic,2), sep = ''),cex = 1.2)
  text(0,-0.6,paste('p = ', signif(r_test$p.value,2), sep = ''),cex = 1.2)
  
  # Inhale onsets
  #text(0,0.55,paste('R = ', round(r_test$statistic,2), sep = ''),cex = 1.2)
  #text(0,0.35,paste('p = ', signif(r_test$p.value,2), sep = ''),cex = 1.2)
  
  # Hartigans' Dip test for unimodality
  #dip.test(angle)
  
  return(data.frame(ID_list,angle,var_angle,r_pval,HR))
}
