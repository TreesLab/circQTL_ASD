#!/usr/bin/env zsh

# FILL IN THESE PATHS
FUSION="<PATH TO FUSION.assoc_test.R>"
LDSC="<PATH TO munge_sumstats.py>"
PLINK="<PATH TO PLINK>"


$LDSC --sumstats daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.gz --daner  --out corrected_daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3
$PLINK --vcf ALL.chr1_GRCh38.genotypes.20170504.vcf.gz --out ALL.GRCh38.genotypes.20170504_SNPs_MAF_chr1 --vcf-half-call m --biallelic-only


Rscript $FUSION \
 --sumstats corrected_daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.sumstats \
 --weights ASD_pcgene_WSLIST_1Mbp_chr1.txt  \
 --weights_dir ./WEIGHTS/chr1/WEIGHTS   \
 --ref_ld_chr ALL.GRCh38.genotypes.20170504_SNPs_MAF_chr \
 --chr 1 \
 --out iPSYCH_ASD_chr1_pcgene_1Mbp.dat




