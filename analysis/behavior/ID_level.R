function (db) {
# input "db" <- output by "block_analyze.R"
  
# Experiment:       respirationCA
# Author:           Martin Grund
# Last update:      November 5, 2019

# FILTER SETTINGS ---------------------------------------------------------

# Number of blocks 
block_min <- 2


# AVERAGE BLOCKS PER PARTICIPANT ------------------------------------------

# Loop participants
for (i in unique(db$ID)) {
    
    db_sub <- subset(db, ID==i & block_filter)
    
    block_num <- nrow(db_sub)
    
    if (block_num >= block_min) {
      
      valid_blocks <- as.integer(paste(sort(db_sub$block), collapse = ''))
      
      results <- c(ID=i, colMeans(db_sub[,c(3:(ncol(db_sub)-6),ncol(db_sub)-1,ncol(db_sub))]), block_num=block_num, valid_blocks=valid_blocks)
    
      if ( tryCatch(is.data.frame(dID1), error=function(cond) FALSE) ) {
        
        dID1[nrow(dID1)+1,] <- results
        
      } else {
        dID1 <- data.frame(t(results))
      }
      
    }
}

return(dID1)

}
