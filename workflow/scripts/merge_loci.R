# Merge SNPs into loci 
## A region of ±500 kb was defined around each variant with a stage I P value below 1 × 10−5. 
## These regions were then merged using bedtools to define nonoverlapping regions.

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)

## Snakemake
### Input
path_snps <- snakemake@input[["snps"]]

### Output
path_unmerged <- snakemake@output[["unmerged_bed"]]
outfile <- snakemake@output[["outfile"]]

## Make .bed file
snp_list_out <- read_tsv(path_snps)

snp_list_out %>%
  mutate(
    chrom = paste0("chr", CHR),
    start = BP - 500000,
    end = BP + 500000,
    start = ifelse(start < 0, 1, start),
  ) %>%
  select(chrom, start, end, SNP, STUDY) %>%
  unite(snp, SNP, STUDY) %>%
  write_tsv(path_unmerged, col_names = FALSE)

## Merge loci
region_locus <- glue("bedtools merge -i {path_unmerged}") %>%
  system(intern = TRUE) %>%
  I() %>%
  read_tsv(col_names = FALSE) %>%
  rename(chr = X1, start = X2, end = X3) %>%
  mutate(locus = glue("{chr}:{start}-{end}"),
         chr = str_replace(chr, "chr", ""),
         chr = as.numeric(chr))

write_csv(region_locus, outfile)