### Get cytogenic bands

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
library(httr)

## Snakemake
### Output
outfile = snakemake@output[['out']]

bands <-
  httr::POST(url = "https://genome.ucsc.edu/cgi-bin/hgTables",
             body = list(clade = "mammal", org = "Human", db = "hg19",
                         hgta_group = "map", hgta_track = "cytoBand",
                         hgta_table = "cytoBand", hgta_regionType = "genome",
                         hgta_outputType = "primaryTable",
                         hgta_doTopSubmit = "get+output")) %>%
  httr::content() %>%
  I() %>%
  read_tsv(col_types = "ciicc") %>%
  dplyr::rename(chr = "#chrom", band = "name") %>%
  dplyr::filter(str_detect(chr, "chr[0-9MTXY]+$")) %>%
  dplyr::mutate(chr = gsub("chr", "", chr),
                chr = gsub("X|chrX", "23", chr),
                chr = gsub("Y|chrY", "24", chr),
                chr = gsub("MT|chrMT|M|chrM", "25", chr),
                chr = as.numeric(chr)) %>%
  tidyr::unite(cytoband, chr, band, sep = "", remove = FALSE)

write_csv(bands, outfile)