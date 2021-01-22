rm(list=ls())
setwd("Directory of the files")

### input format 
# Covariates sample_1 sample_2 ... sample_N
Cov <- read.table("Covariate.txt",sep='\t',header=TRUE,row.names=1)

### input format 
# ID sample_1 sample_2 ... sample_N
snp <- read.table("SNP.M.txt",sep='\t',header=FALSE)
circ <- read.table("circ.M.txt",sep='\t',header=FALSE)
gene <- read.table("gene.M.txt",sep='\t',header=FALSE)

### Read information of covariates
Region <- as.numeric(Cov[1,])
Diagnosis <- as.numeric(Cov[2,])
Age <- as.numeric(Cov[3,])
Sex <- as.numeric(Cov[4,])
RIN <- as.numeric(Cov[5,])
SeqBatch <- as.numeric(Cov[6,])
BrainBank <-  as.numeric(Cov[7,])
PMI <- as.numeric(Cov[8,])

dat <- matrix(0,nrow=nrow(gene),11)
library(mediation)

### Run Mediation tests
for (i in 1:nrow(gene)){
  snp1 <- as.numeric(snp[i,2:106])
  circ1 <- as.numeric(circ[i,2:106])
  gene1 <- as.numeric(gene[i,2:106])
  
  if ((sum(snp1)!=0)&&(sum(gene1!=0))) { 
    lmodel1 <- lm(circ1 ~ snp1 + Region + Diagnosis + Age + Sex + RIN + PMI + SeqBatch + BrainBank)
    lmodel2 <- lm(gene1 ~ snp1 + circ1 + Region + Diagnosis + Age + Sex + RIN + PMI + SeqBatch + BrainBank)
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
