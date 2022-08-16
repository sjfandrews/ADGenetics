library(tidyverse)
library(glue)
`%nin%` = negate(`%in%`)

setwd("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/")

ad_loci.raw <- read_csv("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/ad_loci.csv")

ad_loci <- ad_loci.raw %>%
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
  distinct(locus, .keep_all = T) 
