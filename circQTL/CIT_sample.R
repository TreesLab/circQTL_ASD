rm(list=ls())
setwd("Directory of the files")

### input format 
# ID sample_1 sample_2 ... sample_N
snp <- read.table("SNP.C.txt",sep='\t',header=FALSE)
circ <- read.table("gene.C.txt",sep='\t',header=FALSE)
gene <- read.table("Trait.C.txt",sep='\t',header=FALSE)

dat <- matrix(0,nrow=nrow(gene),8)
library(cit)

### Run Causal inference tests (CIT)
for (i in 1:nrow(gene)){
  snp1 <- as.numeric(snp[i,2:ncol(snp)])
  circ1 <- as.numeric(circ[i,2:ncol(snp)])
  gene1 <- as.numeric(gene[i,2:ncol(snp)])
  if ((sum(snp1)!=0)&&(sum(gene1!=0))) { 
    
    result = cit.bp(snp1,circ1,gene1)
    
    dat[i,1] <- as.vector(snp[i,1])
    dat[i,2] <- as.vector(circ[i,1])
    dat[i,3] <- as.vector(gene[i,1])
    dat[i,4] <- as.numeric(result[1]) # p_cit
    dat[i,5] <- as.numeric(result[2]) # p_TassocL
    dat[i,6] <- as.numeric(result[3]) # p_TassocGgvnL
    dat[i,7] <- as.numeric(result[4]) # p_GassocLgvnT
    dat[i,8] <- as.numeric(result[5]) # p_LindTgvnG
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
  }  
  
  if (i%%1000==0)
  {print(i)}
  gc()
}


### Output
write.csv(dat, file = "out.CIT.csv", row.names = T)

############################################################################################

# Sample Size
ss = ncol(snp)-1

set.seed(1)
n.perm = 100
perm.index = matrix(NA, nrow=ss, ncol=n.perm )
for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:ss )

n.tests_trans <- nrow(gene)
myresults_cit_trans = vector('list', n.tests_trans)
myresults_cit_fdr <- c()

dat <- matrix(0,nrow=nrow(gene),23)

### Run Causal inference tests (CIT) with permutation (n.perm)
library(cit)
for (i in 1:nrow(gene)){
  snp1 <- as.numeric(snp[i,2:ncol(snp)])
  circ1 <- as.numeric(circ[i,2:ncol(snp)])
  gene1 <- as.numeric(gene[i,2:ncol(snp)])
  
  myresults_cit_trans[[i]] = cit.bp(snp1,circ1,gene1, perm.index = perm.index, n.perm = n.perm)
  fdr_cit_one <- fdr.cit(myresults_cit_trans[i])
  
  dat[i,1] <- as.vector(snp[i,1])
  dat[i,2] <- as.vector(circ[i,1])
  dat[i,3] <- as.vector(gene[i,1])
  dat[i,4] <- fdr_cit_one$p.cit
  dat[i,5] <- fdr_cit_one$q.cit
  dat[i,6] <- fdr_cit_one$q.cit.ll
  dat[i,7] <- fdr_cit_one$q.cit.ul
  dat[i,8] <- fdr_cit_one$q.TaL
  dat[i,9] <- fdr_cit_one$q.ll.TaL
  dat[i,10] <- fdr_cit_one$q.ul.TaL
  dat[i,11] <- fdr_cit_one$q.TaGgvL
  dat[i,12] <- fdr_cit_one$q.ll.TaGgvL
  dat[i,13] <- fdr_cit_one$q.ul.TaGgvL
  dat[i,14] <- fdr_cit_one$q.GaLgvT
  dat[i,15] <- fdr_cit_one$q.ll.GaLgvT
  dat[i,16] <- fdr_cit_one$q.ul.GaLgvT
  dat[i,17] <- fdr_cit_one$q.LiTgvG
  dat[i,18] <- fdr_cit_one$q.ll.LiTgvG
  dat[i,19] <- fdr_cit_one$q.ul.LiTgvG
  dat[i,20] <- fdr_cit_one$p_TassocL
  dat[i,21] <- fdr_cit_one$p_TassocGgvnL
  dat[i,22] <- fdr_cit_one$p_GassocLgvnT
  dat[i,23] <- fdr_cit_one$p_LindTgvnG
  
  if (i%%1000==0) 
  {print(i)}
}

### Output
write.csv(dat, file = "out.CIT.perm.csv", row.names = T)