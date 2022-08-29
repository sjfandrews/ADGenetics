# Combine annotation datasets

library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)

# Snakemake 
## Input
# adgwas.path = snakemake@input[['gwas_meta']]

adgwas.path = "intermediate/adgwas_snps.csv"
other.path = "intermediate/other_loci.csv"
tags.path = "intermediate/adgwas_tags.csv"
locus.path = "intermediate/adgwas_loci.csv"
bands.path = "intermediate/cytobands.csv" 
dbsnp.path =  "intermediate/adgwas_dbsnp.csv"
gnomad.path = "intermediate/adgwas_gnomad.csv"
gencode.path = "intermediate/adgwas_gencode.csv"
gvc.path = "intermediate/adgwas_gvc.csv"
snpeff.path = "intermediate/adgwas_snpeff.csv"

## Output
outfile = snakemake@output[['outfile']]

# Import datasets
message(
  "\nImporting...", 
  "\n\t", adgwas.path,
  "\n\t", other.path,
  "\n\t", tags.path,
  "\n\t", locus.path,
  "\n\t", bands.path,
  "\n\t", dbsnp.path,
  "\n\t", gnomad.path,
  "\n\t", gencode.path,
  "\n\t", gvc.path, 
  "\n"
        )

other <- read_csv(other.path, col_types = list(A2 = col_character()))
snps <- read_csv(adgwas.path)
tags <- read_csv(tags.path)
locus <- read_csv(locus.path)
bands <- read_csv(bands.path)
dbsnp <- read_csv(dbsnp.path)
gnomad <- read_csv(gnomad.path)
gencode <- read_csv(gencode.path)
gvc <- read_csv(gvc.path)
snpeff <- read_csv(snpeff.path)

## Merge regions, LD, and cytogenic bands
loci <- tags %>%
  mutate(snp_start = POS, snp_end = POS + 1) %>%
  fuzzyjoin::genome_left_join(locus,
    by = c("CHR" = "chr", "snp_start" = "start", "snp_end" = "end")) %>%
  fuzzyjoin::genome_left_join(bands,
    by = c("CHR" = "chr", "snp_start" = "chromStart",
           "snp_end" = "chromEnd")) %>%
  select(SNP = RS_Number, POS, Alleles, Details,
         tag, tag_clump, locus_ld, locus, cytoband)

## Join datasets together
out <- other %>%
  filter(str_detect(GENE, "APOE")) %>%
  bind_rows(snps, .) %>%
  left_join(select(dbsnp, SNP, dbsnp_gene, Major, Minor, MAF, AF),
            by = "SNP") %>%
  left_join(select(gnomad, SNP, gnomad_minor, gnomad_af, gnomad_maf),
            by = "SNP") %>%
  dplyr::rename(global_maf = MAF, global_af = AF) %>%
  mutate(
    global_maf = ifelse(SNP == "rs139643391", 0.087446, global_maf),
    Major = ifelse(str_detect(SNP, "APOE"), A1, Major),
    Minor = ifelse(str_detect(SNP, "APOE"), A2, Minor),
    FRQ = ifelse(Minor != A2 & nchar(Minor) == 1 & !between(FRQ, 0.42, 0.58),
                 1 - FRQ, FRQ),
    OR = ifelse(Minor != A2 & nchar(Minor) == 1  & !between(FRQ, 0.42, 0.58),
                1 / OR, OR),
    A1 = ifelse(Minor != A2 & nchar(Minor) == 1  & !between(FRQ, 0.42, 0.58),
                A2, A1),
    A2 = ifelse(Minor != A2 & nchar(Minor) == 1  & !between(FRQ, 0.42, 0.58),
                Minor, A2)) %>%
  left_join(select(loci, SNP, locus, cytoband, locus_ld), by = c("SNP")) %>%
  relocate(locus, cytoband, locus_ld) %>%
  left_join(select(gencode, SNP, gencode_gene, gencode_dist, direction),
            by = "SNP") %>%
  left_join(select(gvc, SNP, gvc_gene, gvc_dist), by = "SNP") %>%
  left_join(snpeff, by = c("SNP" = "snp")) %>%
  mutate(locus_ld = as.numeric(locus_ld)) #%>%
  #mutate(locus_ld = ifelse(locus_ld >= 88, locus_ld + 1, locus_ld))

# Export
message("\nExporting....", outfile, "\n")
write_csv(out, outfile)
