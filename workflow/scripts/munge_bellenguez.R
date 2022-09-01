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

data(sumstatsColHeaders)

# Snakemake input & output
bellenguez.path <- snakemake@input[['bellenguez']]
bellenguez_clean.path = snakemake@output[['bellenguez_clean']]
bellenguez_munge.path = snakemake@output[['bellenguez_munge']]

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


### APOE locus added manually from Kunkle summary statistics as not reported

message("\nExtracting Bellenguez 2022 .pdf\n")
bellenguez_raw <- extract_tables(bellenguez.path)

message("\nMunging Bellenguez 2022\n")
bellenguez_clean <- bind_rows(
  bellenguez_raw[[1]] %>% 
    as_tibble() %>%
    slice(4:nrow(.)) %>%
    separate(col = "V1", into = c("snp", "chr", "pos", "gene", 'locus', 'alleles', 'maf', 'or', 'ci', 'p'), 
             sep = " "), 
  
  bellenguez_raw[[2]] %>% 
    as_tibble() %>%
    slice(6:nrow(.)) %>%
    mutate(V1 = str_replace(V1, "IGH gene cluster", "IGH_gene_cluster")) %>%
    separate(col = "V1", into = c('locus', "snp", "chr", "pos", "gene", 'alleles', 'maf', 'or', 'ci', 'p'), 
             sep = " ")
) %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("effectAllele", "otherAllele"), sep = "/") %>%
  mutate(p = str_replace_all(p, "!×!10−| ×!10−|x 10−|!× !10−|!× 10−| × 10−", "e-"),
         or = as.numeric(or), 
         l95 = as.numeric(l95), 
         u95 = as.numeric(u95), 
         p = as.numeric(p), 
         chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         beta = log(or), 
         se = calc_se_from_beta_and_l95(or, l95)) %>%
  add_row(snp = "rs429358", chr = 19, pos = 44908684, gene = "APOE", locus = "APOE", 
          effectAllele = "C", otherAllele = "T", maf = '0.216', 
          or = 3.32, l95 = 3.20, u95 = 3.45, p = 1.2e-881, beta = 1.199, se = 0.019) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  arrange(CHR, BP)

message("\nStandardizing & Harmonizing Bellenguez 2022\n")
bellenguez_munged <- MungeSumstats::format_sumstats(path=bellenguez_clean, 
                                                    ref_genome="GRCh38", 
                                                    convert_ref_genome = "GRCh37",
                                                    allele_flip_check = FALSE, 
                                                    allele_flip_drop = FALSE, 
                                                    bi_allelic_filter = FALSE,
                                                    return_data = TRUE, 
                                                    log_folder = "data/MungeSumstats", 
                                                    force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

message("\nJoining DataFrames\n")
### update BP for SNPs not lifted over  
bellenguez <- bellenguez_munged %>% 
  bind_rows(anti_join(bellenguez_clean, ., by = "SNP")) %>%
  mutate(BP = case_when(
    SNP == "rs139643391" ~ 203743439,
    SNP == "rs6605556" ~ 32583099,
    SNP == "rs1160871" ~ 28168746,
    SNP == "rs35048651" ~ 1631341,
    SNP == "rs199515" ~ 44856641,
    SNP == "rs149080927" ~ 1854258,
    SNP == "rs587709" ~ 54771451,
    TRUE ~ BP
  )) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

print(bellenguez, n = Inf)

## Export Cleaned datasets 
message("\nExporting: ", bellenguez_munge.path, " & ", 
        bellenguez_clean.path, "\n")

write_csv(bellenguez_clean, bellenguez_clean.path)
write_csv(bellenguez, bellenguez_munge.path)

































































