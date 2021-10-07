rm(list=ls())
setwd("Directory of the files")

### input format 
# ID sample_1 sample_2 ... sample_N
snp <- read.table("SNP.M.txt",sep='\t',header=FALSE)
circ <- read.table("circ.M.txt",sep='\t',header=FALSE)
gene <- read.table("gene.M.txt",sep='\t',header=FALSE)

dat <- matrix(0,nrow=nrow(gene),11)
library(mediation)

### Run Mediation tests
for (i in 1:nrow(gene)){
  print(i)
  snp1 <- as.numeric(snp[i,2:ncol(snp)])
  circ1 <- as.numeric(circ[i,2:ncol(snp)])
  gene1 <- as.numeric(gene[i,2:ncol(snp)])
  
  if ((sum(snp1)!=0)&&(sum(gene1!=0))) { 
    lmodel1 <- lm(circ1 ~ snp1 )
    lmodel2 <- lm(gene1 ~ snp1 + circ1)
    result <- mediate(lmodel1, lmodel2, sims=1000, treat = "snp1", mediator = "circ1")
    
    dat[i,1] <- as.vector(snp[i,1])
    dat[i,2] <- as.vector(circ[i,1])
    dat[i,3] <- as.vector(gene[i,1])
    dat[i,4] <- summary(result)$d0 #ACME Beta
    dat[i,5] <- summary(result)$d0.p #ACME pvalue
    dat[i,6] <- summary(result)$z0 #ADE Beta
    dat[i,7] <- summary(result)$z0.p #ADE pvalue
    dat[i,8] <- summary(result)$tau.coef #Total Beta
    dat[i,9] <- summary(result)$tau.p #Total pvalue
    dat[i,10] <- summary(result)$n0 #Prop
    dat[i,11] <- summary(result)$n0.p #Prop pvalue
  }
  else {
    dat[i,1] <- as.vector(snp[i,1])
    dat[i,2] <- as.vector(circ[i,1])
    dat[i,3] <- as.vector(gene[i,1])
    dat[i,4] <- "NA"
    dat[i,5] <- "NA"
    dat[i,6] <- "NA"
    dat[i,7] <- "NA"
    dat[i,8] <- "NA"
    dat[i,9] <- "NA"
    dat[i,10] <- "NA"
    dat[i,11] <- "NA"
  }  
  if (i%%1000==0) 
  {print(i)}
}

### Output
write.csv(dat, file = "out.mediation.csv", row.names = T)