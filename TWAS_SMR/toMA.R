library(tidyverse)
library(readr)


ASD_GWAS = read_delim("daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.gz", delim="\t")

ASD_GWAS.1 = ASD_GWAS %>% 
  dplyr::select(SNP, A1, A2, FRQ_A_18381, OR, SE, P) %>% 
  mutate("b"=log(OR)) %>% 
  mutate("n"= 18381+27969) %>% 
  dplyr::rename("freq"=FRQ_A_18381, "se"=SE, "p"=P) %>%
  dplyr::select(SNP, A1, A2, freq, b, se, p, n) %>% 
  filter(str_detect(SNP, "rs"))
write_delim(ASD_GWAS.1, "daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.ma", delim = "\t")

