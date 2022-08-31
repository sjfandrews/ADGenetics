# GWAS Meta-data 

library(tibble)
library(dplyr)
library(readr)

## Snakemake input/output
outfile <- snakemake@output[["out"]]

## Effective sample size
## https://github.com/neurogenomics/MungeSumstats/blob/master/R/compute_sample_size_neff.R

adgwas <- tribble(
  ~study, ~year, ~neff, ~n, ~n_pca, ~n_ca, ~n_cn, ~n_loci,
    ~ancestry, ~cohorts, ~notes,
  "Lambert", 2013, NA, 74046, NA, 17008, 37154, 19, NA, NA, NA,
  "Kunkle", 2019, NA, 94437, NA, 35274, 59163, 24, NA, NA, NA,
  "Marioni", 2018, NA, 377012, NA, (27696 + 14338 + 25580), (37154 + 245941),
    26, NA, NA, NA,
  "Jansen", 2019, NA, 455258, NA, (24087 + 47793), (55058 + 328320),
    29, NA, NA, NA,
  "Wightman", 2021, NA, 1126563, NA, 90338, 1036225, 38, NA, NA, NA,
  "Bellenguez", 2022, NA, 788989, 46828, (39106 + 46828 + 25392),
    (401577 + 276086), 75, NA, NA, NA,
  "Naj", 2021, NA, NA, NA, (25170 + 20301 + 35214), (41052 + 21839 + 180791),
    29, NA, NA, "29 in stage 1 + 2; 24 w/ UKB in stage 3",
  "Corder", 1993, 0, 1, NA, 1, 1, NA, NA, NA, NA,
  "Jonsson", 2012, 0, 1, NA, 1, 1, NA, NA, NA, NA,
) %>%
  mutate(neff = 4 / ((1 / n_ca) + (1 / n_cn)))

write_csv(adgwas, outfile)
