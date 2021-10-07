### circQTL.precojo.1.txt format:
### SNP_ID A1 A2 freq Beta SE p-value Sample_N

Directory/gcta_1.93.2beta/gcta64 --bfile Directory/plink_files_chr${chr} --chr ${chr} --cojo-file circQTL.precojo.1.txt --cojo-slct --cojo-p 0.05 --out circQTL.cojo.1.txt
