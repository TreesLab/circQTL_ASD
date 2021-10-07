#### 1. Identification of eQTL/circQTL from expression profiles (R code)

One R script: QTL_sample.R

Input:

    (1) Genotyping profiles of SNPs for evaluation of QTL (SNP.txt) 
    (2) The position of each SNP from (2) (snpsloc.txt) 
    (3) Gene expression profiles of genes fro evaluation of QTL (mRNA_expression.txt)
    (4) The position of each gene from (4) (geneloc.txt)

Output:

    out_eqtl.csv

#### 2. Identification of eQTL/circQTL from expression profiles by permutation (R code)

One R script: QTL_permutation_sample.R

Input:

    (1) Genotyping profiles of SNPs for evaluation of QTL (circQTL.1.txt) 
    (2) Gene expression profiles of genes for evaluation of QTL (circRNA.1.txt)

Output:

    circQTL.empP.1.csv

#### 3. Identification of independent circQTL by GCTA-COJO (Shell script)

One shell script: cojo.sh

Input:

    (1) Input file of SNP for cojo (circQTL.precojo.1.txt) 
    (2) plink bfiles

Output:

    circQTL.cojo.1.txt

#### 4. Examination of mediation effects of SNP-circRNA-gene axes (R code)

One R script: Mediation_sample.R

Input:

    (1) Genotyping profiles of SNPs for evaluation of mediation (SNP.M.txt) 
    (2) Expression profiles of circRNAs for evaluation of mediation (circ.M.txt) 
    (3) Expression profiles of genes for evaluation of mediation (gene.M.txt)

Output:

    out.mediation.csv

#### 5. Examination of causal inference effects of SNP-gene-Trait axes (R code)

One R script: CIT_sample.R

Input:

    (1) Genotyping profiles of SNPs for evaluation of CIT (SNP.C.txt) 
    (2) Expression profiles of genes for evaluation of CIT (gene.C.txt) 
    (3) Trait information of each individual (Trait.C.txt)

Output:

    out.CIT.csv
    out.CIT.perm.csv
