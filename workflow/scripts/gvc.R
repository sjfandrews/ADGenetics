# AD Gene verification committee
# https://adsp.niagads.org/index.php/gvc-top-hits-list/

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(googlesheets4)
library(readr)
library(forcats)

## Snakemake
### Input
path_snps <- snakemake@input[["snps"]]
path_gencode <- snakemake@input[["gencode"]]

### Output
outfile <- snakemake@output[["outfile"]]

## Functions
gene_dist <- Vectorize(function(pos, start, end) {
  if (between(pos, start, end)) return(0)
  rg <- c(start, end)
  min(abs(rg - pos))
})

nearest_gene_df <- function(chrom, position, genes) {
  genes %>%
    filter(chr == chrom) %>%
    mutate(dist = gene_dist(position, start, end)) %>%
    arrange(dist) %>%
    slice_head %>%
    mutate(
      direction = ifelse(
        dist == 0,
        "intragenic",
        ifelse(position > end, "downstream", "upstream")),
      pos = position) %>%
    select(-start, -end)
}

## Reading in
message("\nImporting...\n\t", path_snps, "\n\t", path_gencode, "\n")
snps <- read_csv(path_snps)

## Gencode
gencode <- rtracklayer::import(path_gencode) %>%
  as_tibble %>%
  filter(gene_type == "protein_coding" & type == "gene" &
         source == "HAVANA" & !is.na(hgnc_id)) %>%
  select(chr = seqnames, start, end, gene_name) %>%
  mutate(chr = as.character(chr)) %>%
  mutate(chr = gsub("chr", "", chr),
         chr = gsub("X|chrX", "23", chr),
         chr = gsub("Y|chrY", "24", chr),
         chr = gsub("MT|chrMT|M|chrM", "25", chr),
         chr = as.numeric(chr)
  ) %>%
  filter(!is.na(chr))

message("\nImporting...\n\t", "GVC\n")
gvc_url <- paste("https://docs.google.com",
  "spreadsheets/d/166P7ThONaIDyPh2luJHoXDxlnhpsr56a8c_OE_uPQ3w",
  "edit?usp=sharing", sep = "/")

### AD loci
gvc_loci <- read_sheet(gvc_url, sheet = 1, skip = 1) %>%
  janitor::clean_names()

### Candidate genes loci
nearest_gvc_df <- function(chrom, position, genes) {
  genes %>%
    # filter(chr == chrom) %>%
    mutate(dist = gene_dist(position, start, end)) %>%
    arrange(dist) %>%
    slice_head %>%
    mutate(
      direction = ifelse(
        dist == 0,
        "intragenic",
        ifelse(position > end, "downstream", "upstream")),
      pos = position) %>%
    select(-start, -end)
}

gvc_genes_raw <- read_sheet(gvc_url, sheet = 2, skip = 1) %>%
  janitor::clean_names()

gvc_genes <- gvc_genes_raw %>%
  janitor::clean_names() %>%
  mutate(location_hg38 = str_replace_all(location_hg38, ",|chr", "")) %>%
  separate(location_hg38, c("chr38", "start38", "end38")) %>%
  mutate(chr38 = as.numeric(chr38),
         start38 = as.numeric(start38)) %>%
  arrange(chr38, start38) %>%
  left_join(gencode, by = c("gene" = "gene_name")) %>%
  rename(gene_name = gene) %>%
  select(chr, start, end, gene_name)

loci_gvc <- snps %>%
  rename(chr = CHR, pos = BP) %>%
  distinct(SNP, chr, pos) %>%
  left_join(
    purrr::map2_dfr(.$chr, .$pos, nearest_gvc_df, gvc_genes),
    by = c("chr", "pos")) %>%
  mutate(
    gene_name = ifelse(dist > 500000, NA, gene_name),
    dist = ifelse(dist > 500000, NA, dist)) %>%
  rename(gvc_gene = gene_name, gvc_dist = dist)

write_csv(loci_gvc, outfile)