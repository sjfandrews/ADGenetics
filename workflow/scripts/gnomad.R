# Annotate with gnomAD minor allele frequency

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
# library(vcfR)

## Snakemake
### Input
path_gnomad <- snakemake@input[["gnomad"]]
path_snps <- snakemake@input[["snps"]]

### Output
outfile <- snakemake@output[["outfile"]]

## SNPs
snp_list_out <- read_tsv(path_snps)

## gnomAD allele frequencies
gnomadaf_raw <- vcfR::read.vcfR(path_gnomad)

gnomad_af <- gnomadaf_raw %>%
  vcfR::vcfR2tidy(., info_only = TRUE) %>%
  magrittr::extract2("fix") %>%
  select(CHROM, POS, ID, REF, ALT, AF, AF_raw,
         AF_afr, AF_eas, AF_amr, AF_nfe, AF_asj, AF_fin) %>%
  mutate(AF = as.numeric(AF),
         gnomad_maf = ifelse(AF > 0.5, 1 - AF, AF),
         gnomad_minor = ifelse(AF > 0.5, REF, ALT),
         CHROM = as.numeric(CHROM)) %>%
  janitor::clean_names()

gnomad_maf_raw <- gnomad_af %>%
  select(chrom, pos, id, ref, alt, gnomad_minor, af, gnomad_maf) %>%
  add_row(chrom = 7, pos = 28168746, id = "rs1160871", ref = "GTCTT",
          alt = "G", af = 0.4076, gnomad_maf = 0.4076)

gnomad_maf <- bind_rows(
  snp_list_out %>%
    left_join(gnomad_maf_raw,
              by = c("CHR" = "chrom", "BP" = "pos",
                     "A1" = "ref", "A2" = "alt")) %>%
    filter(!is.na(id)),

  snp_list_out %>%
    filter(SNP %in% c("rs9271058", "rs9271192",
                      "rs35048651", "rs540800940")) %>%
    left_join(select(gnomad_maf_raw, id, af, gnomad_maf), by = c("SNP" = "id"))
) %>%
  arrange(CHR, BP) %>%
  select(SNP, CHR, BP, A1, A2, gnomad_minor, gnomad_af = af, gnomad_maf)

write_csv(gnomad_maf, outfile)