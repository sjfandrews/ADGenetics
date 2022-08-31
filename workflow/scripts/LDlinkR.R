# LD Structure
## Assign variants to index snps

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
# library(LDlinkR)

## Snakemake
### Input
path_snps <- snakemake@input[["snps"]]

### Output
out_ld <- snakemake@output[["ld"]]
out_tags <- snakemake@output[["tags"]]

## Functions
### Calculated LD from LDlink
calc_ld <- function(snp_df) {
  message("Computing ld for chr ", snp_df %>% slice(1) %>% pull(chr), "...")

  if (nrow(snp_df) > 1) {
    out <- LDlinkR::LDmatrix(snps = snp_df$SNP,
                             pop = "EUR", r2d = "r2",
                             token = "dbae3c3dc2d4"
    )
  } else {
    out <- NULL
  }
  out
}

### Taging SNPs from LDlink
tag_snps <- function(snp_df) {
  message("Computing ld for chr ", snp_df %>% slice(1) %>% pull(chr), "...")

  if (nrow(snp_df) > 1) {
    out <- LDlinkR::SNPclip(snps = snp_df$SNP,
                            pop = "EUR",
                            r2_threshold = "0.1",
                            maf_threshold = "0.001",
                            token = "dbae3c3dc2d4"
    )
  } else {
    out <- NULL
  }
  out
}

### LD Structure
snps <- read_csv(path_snps)

ld <- snps %>%
  select(SNP, CHR, BP, A1, A2, study) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  mutate(chr = CHR) %>%
  nest(data = c(SNP, chr, BP, A1, A2, study)) %>%
  arrange(CHR) %>%
  mutate(ld = map(data, calc_ld),
         tags = map(data, tag_snps))

write_rds(ld, out_ld)

## Taging SNPs
tags_prune <- ld %>%
  select(tags) %>%
  unnest(cols = tags) %>%
  mutate(tag = str_extract(Details, "rs[:digit:]+"),
         tag = ifelse(is.na(tag), RS_Number, tag)) %>%
  separate(Position, c("CHR", "POS")) %>%
  full_join(ld %>% select(data) %>% unnest(data) %>% select(SNP, chr, BP),
            by = c("RS_Number" = "SNP")) %>%
  mutate(CHR = str_replace(CHR, "chr", ""),
         CHR = as.numeric(CHR),
         POS = as.numeric(POS),
         # Fill in infromation from variants not in 1KG
         CHR = ifelse(is.na(CHR), chr, CHR),
         POS = ifelse(is.na(POS), BP, POS),
         tag = ifelse(is.na(tag), RS_Number, tag),
         tag = as_factor(tag)) %>%
  select(-chr, -BP) %>%
  arrange(CHR, POS) %>%
  mutate(
    locus_ld = lvls_revalue(tag, as.character(seq_along(fct_unique(tag))))
  )

# Clumping
tags <- snps %>%
  # get the minimum P value for each variant
  arrange(CHR, BP, P) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  select(SNP, P) %>%
  # Join the minimum P value for each variant onto the tag df
  left_join(tags_prune, ., by = c("RS_Number" = "SNP")) %>%
  # Find the minimum P value for each ld locus
  group_by(locus_ld) %>%
  arrange(P) %>%
  distinct(locus_ld, .keep_all = TRUE) %>%
  ungroup %>%
  # Get the tagging SNP by clump and join back onto the tag df
  select(tag_clump = RS_Number, locus_ld) %>%
  left_join(tags_prune, ., by = "locus_ld")

write_csv(tags, out_tags)