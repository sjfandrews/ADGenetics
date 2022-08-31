# Annotate with Nearest Gene
## Nearest protein-coding gene according to GENCODE release 40

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
# library(rtracklayer)

# Snakemake
## Input
path_snps <- snakemake@input[["snps"]]
path_gencode <- snakemake@input[["gencode"]]

## Output
outfile <- snakemake@output[["outfile"]]

# Functions

gene_dist <- Vectorize(function(pos, start, end) {
  if (between(pos, start, end)) return(0)
  rg <- c(start, end)
  min(abs(rg - pos))
})

nearest_gene_df <- function(chrom, position, genes) {
  message("chr: ", chrom, " pos: ", position)
  out <- genes %>%
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
  out
}

# Gencode v40
message("\nImport: ", path_gencode, "...\n")
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

# Read SNPS
message("\nAnnotate with nearest gene\n")
snps <- read_csv(path_snps)

loci_gene <- snps %>%
  rename(chr = CHR, pos = BP) %>%
  distinct(SNP, chr, pos) %>%
  left_join(
    purrr::map2_dfr(.$chr, .$pos, nearest_gene_df, gencode),
    by = c("chr", "pos")) %>% 
  rename(gencode_gene = gene_name, gencode_dist = dist)

# Export
write_csv(loci_gene, outfile)