
## Kunkle, B. W. et al. Nat Genet 51, 414–430 (2019).
  
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
kunkle.path <- snakemake@input[['kunkle']]

### Output
kunkle_munge.path = snakemake@output[['kunkle_munge']]
kunkle_clean.path = snakemake@output[['kunkle_clean']]

## kunkle 2019
message("\nMunging Kunkle 2019\n")

kunkle_raw <- extract_tables(kunkle.path)
k1 <- as_tibble(kunkle_raw[[1]]) %>% 
  slice(-c(1:7,27,30)) %>%
  select(-V3) %>%
  mutate(V1 = str_replace_all(V1, "!×!10−| ×!10−|x 10−", "e-"), 
         V2 = str_replace_all(V2, "!×!10−| ×!10−|x 10−", "e-")) %>%
  separate(V1, c("snp", 'chr', 'pos', 'gene','alleles', 'maf', 'or1', 'ci1', 'p1', 
                 'or2', 'ci2', 'p2'), sep = " ") %>%
  separate(V2, c('or', 'ci', 'p'), sep = " ") %>% 
  mutate(or = na_if(or, ""),
         or = ifelse(is.na(or), or1, or), 
         ci = ifelse(is.na(ci), ci1, ci), 
         p = ifelse(is.na(p), p1, p)) %>%
  select(snp, chr, pos, gene, alleles, maf, or, ci, p) %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("otherAllele", "effectAllele"), sep = "/")

k2 <- as_tibble(kunkle_raw[[2]]) %>% 
  slice(-c(1:4,14:16)) %>% 
  select(V1, V2, V3, V4, V6, V11, V12) %>%
  mutate(V12 = str_replace_all(V12, "!×!10−| ×!10−|x 10−", "e-")) %>%
  magrittr::set_colnames(c("snp", "chr", "pos", "gene", "maf", "or", "p")) %>%
  separate(gene, c("gene",'alleles'), sep = " ") %>%
  separate(or, c("or",'ci'), sep = " ") %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("otherAllele", "effectAllele"), sep = "/")

kunkle_clean <- bind_rows(
  k1, k2
) %>%
  mutate(p = as.numeric(p),
         or = as.numeric(or), 
         l95 = as.numeric(l95), 
         u95 = as.numeric(u95), 
         beta = log(or), 
         maf = as.numeric(maf), 
         se = calc_se_from_beta_and_l95(or, l95),
         chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         snp = str_replace(snp, "rs7920721g", "rs7920721")) %>%
  arrange(chr, pos) %>% 
  filter(p < 5e-8) %>%
  group_by(snp) %>%
  slice(which.min(p)) %>% 
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

kunkle <- MungeSumstats::format_sumstats(path=kunkle_clean, 
                                          ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                          log_folder = "data/MungeSumstats", 
                                          force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(kunkle_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

## Export Cleaned datasets 
write_csv(kunkle_clean, kunkle_clean.path)
write_csv(kunkle, kunkle_munge.path)
