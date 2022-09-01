
## Lambert, J.-C. et al. Nat Genet 45, 1452–1458 (2013).
  
library(MungeSumstats)
library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(plyr)
library(dplyr)
library(stringr)
library(forcats)
library(tabulizer)
library(fuzzyjoin)

## Misc functions
calc_se_from_beta_and_l95 = function(or, l95){
  logOR = log(or)
  logOR_L95 = log(l95)
  logOR_SE = (logOR-logOR_L95)/1.96
  logOR_SE
}

calc_beta_from_z_n_freq = function(Z, FRQ, N){
  out = Z/sqrt(2*FRQ*(1-FRQ)*(N+Z^2))
  out
}

data(sumstatsColHeaders)

## Snakemake 
### Input
lambert.path <- snakemake@input[['lambert']]

### Output
lambert_munge.path = snakemake@output[['lambert_munge']]
lambert_clean.path = snakemake@output[['lambert_clean']]

## Lambert 2013
message("\nMunging Lambert 2013\n")

lambert_raw <- extract_tables(lambert.path)
lambert_clean <- as_tibble(lambert_raw[[1]]) %>% 
  slice(-c(1,2,3,4,23,34)) %>%
  select(V1:V5, V12:V13) %>%
  magrittr::set_colnames(c("snp", "chr", "pos", "gene", "alleles", "or", "p", "hetpval")) %>%
  separate(alleles, c("alleles", "maf"), sep = " ") %>%
  mutate(snp = na_if(snp, ""), 
         effect = rep(c("or", "ci"), 21)) %>%
  fill(snp, .direction = "down") %>%
  pivot_wider(names_from = effect, values_from = or) %>%
  fill(ci, .direction = "up") %>%
  filter(!is.na(or)) %>%
  mutate(ci = str_replace_all(ci, "\\(|\\)", "")) %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("otherAllele", "effectAllele"), sep = "/") %>%
  mutate(p = str_replace_all(p, " × 10−", "e-"), 
         maf = as.numeric(maf),
         p = as.numeric(p),
         or = as.numeric(or), 
         l95 = as.numeric(l95), 
         u95 = as.numeric(u95), 
         beta = log(or), 
         se = calc_se_from_beta_and_l95(or, l95),
         chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         snp = str_replace(snp, "rs3865444g", "rs3865444"), 
         snp = str_replace(snp, "rs8093731g", "rs8093731")) %>%
  read_sumstats( .,
    standardise_headers = T,
    mapping_file = sumstatsColHeaders
  ) %>%
  filter(P < 5e-8) %>%
  mutate(CHR = as.numeric(CHR)) %>%
  arrange(CHR, BP) %>%
  as_tibble()

lambert <- MungeSumstats::format_sumstats(path=lambert_clean, 
                               ref_genome="GRCh37", 
                               allele_flip_check = FALSE, 
                               allele_flip_drop = FALSE, 
                               bi_allelic_filter = FALSE,
                               return_data = TRUE, 
                               log_folder = "data/MungeSumstats", 
                               force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(lambert_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

## Export Cleaned datasets 
write_csv(lambert_clean, lambert_clean.path)
write_csv(lambert, lambert_munge.path)
