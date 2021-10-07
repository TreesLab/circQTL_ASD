#!/usr/bin/env zsh

# FILL IN THESE PATHS
SMR="<PATH TO smr_Linux>"
PLINK="<PATH TO PLINK>"

$PLINK --vcf ALL.chr1_GRCh38.genotypes.20170504.vcf.gz --out ALL.GRCh38.genotypes.20170504_chr1 --vcf-half-call m --biallelic-only

$SMR \
--bfile ALL.GRCh38.genotypes.20170504_chr1 \
--gwas-summary daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.ma \
--beqtl-summary ASD_pcgene_matrixEQTL \
--out smr_result/ASD_circRNA_chr1 \
--thread-num 4

