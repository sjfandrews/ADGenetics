library(MungeSumstats)
library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(plyr)
library(dplyr)
library(stringr)
library(forcats)
library(data.table)

# Snakemake

## Input
lambert.path <- snakemake@input[['lambert']]
kunkle.path <- snakemake@input[['kunkle']]
marioni.path <- snakemake@input[['marioni']]
jansen.path <- snakemake@input[['jansen']]
wightman.path <- snakemake@input[['wightman']]
bellenguez.path <- snakemake@input[['bellenguez']]
other.path = snakemake@input[['other']]
adgwas.path = snakemake@input[['meta']]

## Output
out_vcf = snakemake@output[['out_vcf']]
out_csv = snakemake@output[['out_csv']]
out_tsv = snakemake@output[['out_tsv']]
out_list = snakemake@output[['out_list']]

# Combine Files 
adgwas <- read_csv(adgwas.path)

other <- read_csv(other.path,
                  col_types = list(A2 = col_character()))

snps <- list(
  Lambert = lambert.path,
  Kunkle = kunkle.path,
  Marioni = marioni.path,
  Jansen = jansen.path,
  Wightman = wightman.path,
  Bellenguez = bellenguez.path) %>%
  map(read_csv) %>%
  map(select, SNP, CHR, BP, A1, A2, any_of(c("GENE")), FRQ, OR, P) %>%
  imap_dfr(~ mutate(.x, study = .y)) %>%
  bind_rows(other %>% filter(!str_detect(GENE, "APOE"))) %>%
  arrange(CHR, BP)

write_csv(snps, out_csv)

## GWAS VCF
snp_list <- snps %>% 
  left_join(adgwas, by = "study") %>% 
  group_by(CHR, BP) %>%
  slice(which.max(neff)) %>%
  ungroup() %>%
  select(SNP, CHR, BP, A1, A2, OR, P, study) %>%
  rename(STUDY = study)

snp_list_munged <- MungeSumstats::format_sumstats(path = snp_list,
                                 ref_genome = "GRCh37",
                                 allele_flip_check = TRUE,
                                 allele_flip_drop = FALSE,
                                 bi_allelic_filter = FALSE,
                                 return_data = TRUE,
                                 log_folder = "data/MungeSumstats",
                                 force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

snp_list_out <- snp_list_munged %>%
  bind_rows(anti_join(snp_list, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) 

snp_list_out %>%
  as.data.table() %>% 
  MungeSumstats::write_sumstats(out_vcf, write_vcf = TRUE)

write_tsv(snp_list_out, out_tsv)

## Unique variants
snp_list <- distinct(snps, SNP)
write_delim(snp_list, out_list, col_names = FALSE)

