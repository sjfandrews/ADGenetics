# Extract MAF from NCBI DBSNP

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
# library(rsnps)

## Snakemake
## Input
path_snps <- snakemake@input[["snps"]]

## Output
outfile <- snakemake@output[["outfile"]]

## Import SNPs
dbsnp_info <- read_delim(path_snps, delim = "\t", col_names = FALSE) %>%
  pull(X1) %>%
  rsnps::ncbi_snp_query()

global_maf <- dbsnp_info %>%
  select(query, gene, maf_population) %>%
  unnest(cols = c(maf_population)) %>%
  filter(study == "dbGaP_PopFreq") %>%
  group_by(query) %>%
  slice(which.max(MAF)) %>%
  ungroup() %>%
  mutate(
    AF = MAF,
    MAF = ifelse(MAF > 0.5, 1 - MAF, MAF),
    Major = ifelse(AF > 0.5, Minor, ref_seq),
    Minor = ifelse(AF < 0.5, Minor, ref_seq),
  ) %>%
  relocate(Major, .before = Minor) %>%
  mutate(gene = na_if(gene, "")) %>%
  rename(dbsnp_gene = gene, SNP = query)

## Export
write_csv(global_maf, outfile)