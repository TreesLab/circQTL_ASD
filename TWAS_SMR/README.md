# Performing transcriptome-wide association studies (TWAS)

# 1. TWAS-FUSION
     (http://gusevlab.org/projects/fusion/)

     Step1: runs_FUSION.compute_weights.sh
           
           Input: 
                  (1) SNP dataset(in PLINK format): ASD.bed; ASD.bim; ASD.fam
                  (2) protein-coding gene expression on chromosome 1 : ASD_pcgene_GEXP_chr1.txt 
                  
           Output:
                  WEIGHTS/         
                  
     Step2: runs_FUSION.assoc_test.sh on chromosome 1
            
            Input: 
                   (1) ASD GWAS summuary statistics (daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.gz),
                       please download from https://www.med.unc.edu/pgc/download-results/
                   
                   (2) A reference sample for LD estimation (1000 Genomes: ALL.chr1_GRCh38.genotypes.20170504.vcf.gz),
                       please download from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
                       
                   (3) ASD_pcgene_WSLIST_1Mbp_chr1.txt    
                       
                   (4) WEIGHTS/
                   
            Output: 
                   iPSYCH_ASD_chr1_pcgene_1Mbp.dat
                                     
                                     
# 2. SMR    
     (https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis)
     
     Step1 : toMA.R
     
            To convert GWAS data into MA format
     
            Input :
                     (1) ASD GWAS summuary statistics (daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.gz),
                       please download from https://www.med.unc.edu/pgc/download-results/
                       
            output:
                     daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.ma          
     
     Step2 : runs_SMR.sh
     
             Input :
                     (1) A reference sample for LD estimation (1000 Genomes: ALL.chr1_GRCh38.genotypes.20170504.vcf.gz),
                         please download from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
             
                     (2) ASD protein coding genes matrixEQTL results (in BESD format): 
                         ASD_pcgene_matrixEQTL.esi; ASD_pcgene_matrixEQTL.epi; ASD_pcgene_matrixEQTL.besd  
              
             Output :
                      ASD_pcgene_chr1.smr  
                         
                                          
