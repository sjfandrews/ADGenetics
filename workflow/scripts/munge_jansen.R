# Jansen 2019 
## Jansen, I. E. et al. Nat Genet 51, 404–413 (2019).
  
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
jansen.path <- snakemake@input[['jansen_table']]
jansen_ss.path <- snakemake@input[['jansen_ss']]

### Output
jansen_clean.path = snakemake@output[['jansen_clean']]
jansen_munge.path = snakemake@output[['jansen_munge']]

## Munging 
message("\nMunging Jansen 2019\n")

jansen_raw <- extract_tables(jansen.path)
jansen_ss <- read_tsv(jansen_ss.path) %>%
  select(-uniqID.a1a2, Nsum, -Neff, -dir, -Z) %>%
  rename(effectAllele = A1, otherAllele = A2)

jansen_clean <- as_tibble(jansen_raw[[1]]) %>%
  slice(-c(1:5, 19)) %>%
  mutate(V1 = str_replace_all(V1, " 2.59 ×", " 2.59E10-4"),
         V1 = str_replace_all(V1, "!×!| ×!|x |!× !|!× | × ", "E"), 
         V1 = str_replace_all(V1, " NA ", " NA NA "), 
         V1 = str_replace_all(V1, "− ", "-")) %>%
  separate(V1, into = c("locus", "chr", "gene", "snp1", "P1", "SNP2", "P2", "snp", "pos", "a2", "a1", "maf", "z", "p", "dir"), sep = " ") %>%
  select(snp, gene) %>%
  left_join(jansen_ss, by = c("snp" = "SNP")) %>%
  rename(n = Nsum) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  filter(P < 5e-8)

jansen <- MungeSumstats::format_sumstats(path=jansen_clean, 
                                         ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(jansen_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE, OR = exp(BETA))

print(jansen, n = Inf)

## Export Cleaned datasets 
write_csv(jansen, jansen_munge.path)
write_csv(jansen_clean, jansen_clean.path)

