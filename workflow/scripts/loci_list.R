library(MungeSumstats)
library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(plyr)
library(dplyr)
library(stringr)
library(forcats)
# library(data.table)

# Snakemake

## Input
# path_gwas <-
#   c("lambert", "kunkle", "marioni", "jansen", "wightman", "bellenguez") %>%
#   sapply(function(x) snakemake@input[[x]])
path_gwas <- list(
  Lambert = snakemake@input[["lambert"]],
  Kunkle = snakemake@input[["kunkle"]],
  Marioni = snakemake@input[["marioni"]],
  Jansen = snakemake@input[[ "jansen"]],
  Wightman = snakemake@input[["wightman"]],
  Bellenguez = snakemake@input[["bellenguez"]])
path_other <- snakemake@input[["other"]]
path_meta <- snakemake@input[["meta"]]

## Output
out_vcf <- snakemake@output[["out_vcf"]]
out_csv <- snakemake@output[["out_csv"]]
out_tsv <- snakemake@output[["out_tsv"]]
out_list <- snakemake@output[["out_list"]]

# Combine Files
meta <- read_csv(path_meta)

other <- read_csv(path_other,
                  col_types = list(A2 = col_character()))

snps <- path_gwas %>%
  map(read_csv) %>%
  map(select, SNP, CHR, BP, A1, A2, any_of(c("GENE")), FRQ, OR, P) %>%
  imap_dfr(~ mutate(.x, study = .y)) %>%
  bind_rows(other %>% filter(!str_detect(GENE, "APOE"))) %>%
  # mutate(study = str_to_title(study)) %>%
  arrange(CHR, BP) 

write_csv(snps, out_csv)

## GWAS VCF
snp_list <- snps %>%
  left_join(meta, by = "study") %>%
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
  data.table::as.data.table() %>%
  MungeSumstats::write_sumstats(., out_vcf, write_vcf = TRUE)

write_tsv(snp_list_out, out_tsv)

## Unique variants
snp_list <- distinct(snps, SNP)
write_delim(snp_list, out_list, col_names = FALSE)
