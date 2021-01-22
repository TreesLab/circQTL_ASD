###Trans-genetic effects of circular RNA expression quantitative trait loci and potential causal mechanisms in autism
---
Te-Lun Mai, Chia-Ying Chen, Tai-Wei Chiang, and Trees-Juen Chuang* (2021) "Trans-genetic effects of circular RNA expression quantitative trait loci and potential causal mechanisms in autism", to be submitted.
---

The scripts used in this study include three parts:

###1. Identification of eQTL/circQTL from expression profiles (R code)

One R script: QTL_sample.R

Input:

    (1) Confounding factors of each individual for all samples (Covariates.txt)
    (2) Genotyping profiles of SNPs for evaluation of QTL (SNP.txt) 
    (3) The position of each SNP from (2) (snpsloc.txt) 
    (4) Gene expression profiles of genes fro evaluation of QTL (mRNA_expression.txt)
    (5) The position of each gene from (4) (geneloc.txt)

Output:

    out_eqtl.csv

###2. Examination of mediation effects of SNP-circRNA-gene axes (R code)

One R script: Mediation_sample.R

Input:

    (1) Confounding factors of each individual for all samples (Covariates.txt)
    (2) Genotyping profiles of SNPs for evaluation of mediation (SNP.M.txt) 
    (3) Expression profiles of circRNAs for evaluation of mediation (circ.M.txt) 
    (4) Expression profiles of genes for evaluation of mediation (gene.M.txt)

Output:

    out.mediation.csv

###3. Examination of causal inference effects of SNP-gene-Trait axes (R code)

One R script: CIT_sample.R

Input:

    (1) Confounding factors of each individual for all samples (Covariates.txt)
    (2) Genotyping profiles of SNPs for evaluation of mediation (SNP.C.txt) 
    (3) Expression profiles of genes for evaluation of mediation (gene.C.txt) 
    (4) Trait information of each individual (Trait.C.txt)

Output:

    out.CIT.csv
    out.CIT.perm.csv
