
## Wightman, D. P. et al. Nat Genet 53, 1276–1282 (2021).
  
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
wightman.path <- snakemake@input[['wightman_table']]
wightman_ss.path <- snakemake@input[['wightman_ss']]

### Output
wightman_clean.path = snakemake@output[['wightman_clean']]
wightman_munge.path = snakemake@output[['wightman_munge']]

## Wightman 2021
message("\nMunging Wightman 2021\n")

wightman_raw <- extract_tables(wightman.path)
wightman_ss <- read_tsv(wightman_ss.path) %>%
  select(-p, -N) %>%
  rename(effectAllele = testedAllele)

wightman_clean <- as_tibble(wightman_raw[[1]]) %>%
  slice(-c(1:3)) %>%
  mutate(V1 = str_replace_all(V1, "!×!10−| ×!10−|x 10−|!× !10−|!× 10−| × 10−", "e-"), 
         V1 = str_replace_all(V1, ",", ""), 
         V2 = str_replace_all(V2, ",", "")) %>%
  separate(V1, into = c("locus", "gene", "pos", "snp", "ea", "freq", "pval"), sep = " ") %>%
  separate(pos, into = c("chr", "pos"), sep = ":") %>%
  rename(n = V2) %>%
  # select(snp, chr, pos, gene, freq, p, n)  %>%
  mutate(chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         freq = as.numeric(freq),
         n = as.numeric(n), 
         pval = as.numeric(pval), 
         pval = ifelse(is.na(pval), 1.0e-300, pval)) %>%
  left_join(wightman_ss, by = c("chr", "pos" = "PosGRCh37")) %>%
  mutate(freq = ifelse(ea != effectAllele, 1-freq, freq)) %>% 
  select(-ea) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) 

wightman <- MungeSumstats::format_sumstats(path=wightman_clean, 
                                         ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(wightman_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = ifelse(is.infinite(Z), 500, Z),
         BETA = calc_beta_from_z_n_freq(Z, FRQ, N), 
         OR = exp(BETA), 
         )

print(wightman, n = Inf)

## Export Cleaned datasets 
write_csv(wightman_clean, wightman_clean.path)
write_csv(wightman, wightman_munge.path)

