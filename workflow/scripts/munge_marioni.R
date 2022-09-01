
## Marioni, R. E. et al. Translational Pscyh. (2018).

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
marioni.path <- snakemake@input[['marioni']]

### Output
marioni_clean.path = snakemake@output[['marioni_clean']]
marioni_munge.path = snakemake@output[['marioni_munge']]

## Marioni 2018
message("\nMunging Marioni 2018\n")

marioni_ss <- read_tsv(marioni.path) %>%
  select(SNP = DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, OR, SE, Z, P) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

marioni <- MungeSumstats::format_sumstats(path=marioni_ss, 
                                         ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(marioni_ss, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) 


## Export Cleaned datasets 
write_csv(marioni, marioni_munge.path)
write_csv(marioni_ss, marioni_clean.path)
