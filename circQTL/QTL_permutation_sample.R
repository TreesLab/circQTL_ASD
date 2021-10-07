rm(list=ls())
setwd("Directory of the files")

library(RcppEigen)

ncirc = 1
nsam = 105

for (i in 1:ncirc){
  
  print(i)
  
  base.dir = getwd();

  snp_file_name = paste(base.dir,"/circQTL.",i,".txt",sep="")
  snp <- read.table(snp_file_name, sep='\t', header=TRUE)
  
  expression_file_name = paste(base.dir,"/circRNA.",i,".txt",sep="") 
  circ <- read.table(expression_file_name, sep='\t', header=TRUE)
  
  ##################################################
  if (nrow(snp)>=1){
    
    dat1 <- matrix(0,1,nrow=nrow(snp))
    for (j in 1:nrow(snp)){
      circ_o <- as.numeric(circ[1,2:(nsam+1)])
      snp_o <- as.numeric(snp[j,2:(nsam+1)])
      lmodel <- lm(circ_o ~ snp_o )
      dat1[j] <- summary( lmodel )$coefficients[2,4] # p-value
    }
    ##################################################
    
    perm = 10000
    dat <- matrix(0,nrow=nrow(snp),perm)

    for (k in 1:perm){
      
      index = sample(nsam,nsam)
      circ_o <- as.numeric(circ[1,2:(nsam+1)])[index]
      

      for (j in 1:nrow(snp)){
        
        snp_o <- as.numeric(snp[j,2:(nsam+1)])
        lmodel <- fastLm(circ_o ~ snp_o )
        
        dat[j,k] <- summary( lmodel )$coefficients[2,4] # p-value
      }
    }
    
    empP <- matrix(0,nrow=nrow(snp),2)
    for (m in 1:nrow(snp)){
      empP[m,1] <- as.vector(snp[m,1])
      empP[m,2] <- (1+sum(dat1[m] > dat[m,]))/(perm+1)
    }
    
    out_file_name = paste(base.dir,"circQTL.empP.",i,".csv",sep="")
    write.csv(empP, file = out_file_name)
  }else{
    out_file_name = paste(base.dir,"circQTL.empP.",i,".csv",sep="")
    write.csv(snp, file = out_file_name)
    
  }
}