# Descriptives of AD GWAS

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
library(glue)
`%nin%` <- negate(`%in%`)

ad_loci <- "/results/ad_loci.csv" %>%
  read_csv() %>%
  filter(study %nin% c("Jonsson", "Reiman")) %>%
  select(locus, cytoband, locus_ld, SNP, CHR, BP, A1, A2, FRQ, OR, P, study)

known_loci <- ad_loci %>%
  filter(study %nin% c("Bellenguez", "Wightman"))

## 42 novel loci discovered by Bellenguez and Whightman
ad_loci %>%
  filter(study %in% c("Bellenguez", "Wightman")) %>%
  anti_join(known_loci, by = "locus") %>%
  # arrange(P) %>%
  # group_by(study) %>%
  distinct(locus, .keep_all = TRUE)


ad_loci %>%
  filter(study %in% c("Bellenguez", "Wightman")) %>%
  group_by(locus)

ad_loci %>%
  filter(study %in% c("Bellenguez", "Wightman")) %>%
  group_by(locus_ld)

ad_loci %>%
  filter(study %in% c("Bellenguez", "Wightman")) %>%
  group_by(study)  %>%
  distinct(locus) %>%
  ungroup() %>%
  count(study)
