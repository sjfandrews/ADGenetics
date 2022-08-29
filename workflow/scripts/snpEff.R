## Extract snpEff annotations 

library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(plyr)
library(dplyr)
library(stringr)
library(forcats)
# library(vcfR)

## Snakemake
### Input
file_ann.path = snakemake@input[['ann']]

### Output
outfile = snakemake@output[['outfile']]

## snpEff
raw_snpeff <- file_ann.path %>%
  vcfR::read.vcfR() %>%
  vcfR::vcfR2tidy(., info_only = FALSE)

message("\nExtract Meta\n")
ann_cols <- raw_snpeff$meta %>%
  filter(Tag == "INFO" & ID == "ANN") %>%
  pull(Description) %>%
  str_match("'(.+)'") %>%
  .[[2]] %>%
  str_split(" \\| ") %>%
  .[[1]]

message("\nExtract gt\n")
ann <- raw_snpeff$gt %>%
  select(ChromKey, POS, SNP = gt_SNP) %>%
  left_join(., raw_snpeff$fix, by = c("ChromKey", "POS")) %>%
  select(-ChromKey, -ID, -QUAL, -FILTER, -LOF, -NMD) %>%
  mutate(ANN = strsplit(ANN, ",")) %>%
  unnest(cols = ANN) %>%
  separate(ANN, into = ann_cols, sep = "\\|") %>%
  janitor::clean_names() %>%
  mutate(chrom = as.numeric(chrom))

message("\nMunge\n")
snpeff <- ann %>%
  filter(transcript_bio_type %in% c("", "protein_coding")) %>%
  group_by(snp) %>%
  slice_head() %>%
  ungroup() %>%
  arrange(chrom, pos)

message("\nExport: ", outfile, "\n")
write_csv(snpeff, outfile)