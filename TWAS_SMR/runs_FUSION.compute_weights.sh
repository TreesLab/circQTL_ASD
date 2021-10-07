#!/usr/bin/env zsh

# FILL IN THESE PATHS
GCTA="<PATH TO GCTA>"
PLINK="<PATH TO PLINK>"
GEMMA="<PATH TO GEMMA>"
FUSION="<PATH TO FUSION.compute_weights.R>"

# Input: (1) SNP 
PRE_GENO="ASD"

# Input : (2) protein-coding gene expression on chromosome 1
PRE_GEXP="ASD_pcgene_GEXP_chr1.txt"

OUT_DIR="./WEIGHTS"

rm -r -f --parents tmp
rm -r -f --parents hsq
rm -r -f --parents out
rm -r -f $OUT_DIR

mkdir --parents tmp
mkdir --parents hsq
mkdir --parents out

# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir $OUT_DIR


cat $PRE_GEXP | head -n 1 | tr '\t' '\n' | tail -n+5  | awk '{print "0" "\t" $0}'> PRE_GEXP.ID 
ve_in=$(cat $PRE_GEXP | wc -l) 


cat $PRE_GEXP | awk -vs=2 -ve=$ve_in 'NR > s && NR <= e' | while read PARM;do                             

 CHR=`echo $PARM | awk  '{print $2}'`
 P0=`echo $PARM | awk '{print $3 - 1e6}'`
 P1=`echo $PARM | awk '{print $4 + 1e6}'`
 GNAME=`echo $PARM| awk '{print $1":chr"$2"_"$3"_"$4}'`

 OUT="tmp/$GNAME"
 
 echo $GNAME $CHR $P0 $P1 >> ASD.pos
 echo $PARM | tr ' ' '\n' | tail -n+5 | paste PRE_GEXP.ID - > $OUT.pheno
 # Get the locus genotypes for all samples and set current gene expression as the phenotype
 $PLINK --bfile $PRE_GENO --pheno $OUT.pheno --make-bed --out $OUT --keep $OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1 --extract $PRE_GENO.bim
 
 # Process all samples together (for reference purposes only since this is mult-ethnic data)
 
 FINAL_OUT="$OUT_DIR/$GNAME"
 
 Rscript $FUSION --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_plink $PLINK --PATH_gcta $GCTA --PATH_gemma $GEMMA --models lasso,top1,enet
 # ALTERNATIVELY ADD COVARIATES HERE USING THE --covar FLAG
 # MINIMAL COMMAND IS: `Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.$pop.tmp --out $FINAL_OUT`


 # Append heritability output to hsq file
 cat $FINAL_OUT.hsq >> hsq/ALL.hsq
 

 # Clean-up just in case
 rm -f $FINAL_OUT.hsq $OUT.tmp.*

done

