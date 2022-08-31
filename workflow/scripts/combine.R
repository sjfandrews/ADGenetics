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

path <- list(
  adgwas = "intermediate/adgwas_snps.csv",
  other = "intermediate/other_loci.csv",
  tags = "intermediate/adgwas_tags.csv",
  locus = "intermediate/adgwas_loci.csv",
  bands = "intermediate/cytobands.csv",
  dbsnp =  "intermediate/adgwas_dbsnp.csv",
  gnomad = "intermediate/adgwas_gnomad.csv",
  gencode = "intermediate/adgwas_gencode.csv",
  gvc = "intermediate/adgwas_gvc.csv",
  snpeff = "intermediate/adgwas_snpeff.csv"
)

if (exists("snakemake")) {
  path <- path %>%
    imap(~ snakemake@input[[.y]])
  outfile <- snakemake@output[["outfile"]]
} else {
  outfile <- "results/adgwas_loci.csv"
}

# Import datasets
read_input <- function(path, input_name) {
  message("\t", path)
  if (input_name == "other") {
    read_csv(path, col_types = list(A2 = col_character()))
  } else {
    read_csv(path)
  }
}

message("Importing...")
input <- imap(path, read_input)
message("")

## Merge regions, LD, and cytogenic bands
loci <- input$tags %>%
  mutate(snp_start = POS, snp_end = POS + 1) %>%
  fuzzyjoin::genome_left_join(input$locus,
    by = c("CHR" = "chr", "snp_start" = "start", "snp_end" = "end")) %>%
  fuzzyjoin::genome_left_join(input$bands,
    by = c("CHR" = "chr", "snp_start" = "chromStart",
           "snp_end" = "chromEnd")) %>%
  select(SNP = RS_Number, POS, Alleles, Details,
         tag, tag_clump, locus_ld, locus, cytoband)

### For updating APOE locus
apoe_locus <- loci %>%
  filter(SNP == "rs429358") %>%
  slice(1) %>%
  select(locus, cytoband, locus_ld)


## Join datasets together
out <- input$other %>%
  filter(str_detect(GENE, "APOE")) %>%
  bind_rows(input$adgwas, .) %>%
  left_join(select(input$dbsnp, SNP, dbsnp_gene, Major, Minor, MAF, AF),
            by = "SNP") %>%
  left_join(select(input$gnomad, SNP, gnomad_minor, gnomad_af, gnomad_maf),
            by = "SNP") %>%
  dplyr::rename(global_maf = MAF, global_af = AF) %>%
  left_join(select(loci, SNP, locus, cytoband, locus_ld), by = c("SNP")) %>%
  relocate(locus, cytoband, locus_ld) %>%
  mutate(
    ## Manual fixing missing data for rs139643391 & rs149080927
    Minor = ifelse(SNP == "rs139643391", gnomad_minor, Minor),
    gnomad_maf = ifelse(SNP == "rs149080927", 0.3829, gnomad_maf),
    ## Updating Reiman APOE loci
    gnomad_maf = ifelse(str_detect(SNP, "APOE"), FRQ, gnomad_maf),
    locus = ifelse(str_detect(SNP, "APOE"), apoe_locus$locus, locus),
    locus_ld = ifelse(str_detect(SNP, "APOE"), apoe_locus$locus_ld, locus_ld),
    cytoband = ifelse(str_detect(SNP, "APOE"), apoe_locus$locus_ld, cytoband),
    Major = ifelse(str_detect(SNP, "APOE"), A1, Major),
    Minor = ifelse(str_detect(SNP, "APOE"), A2, Minor),
  ) %>%
  mutate(
    ## Checking Allele filps for frequncy
    FRQ = ifelse(Minor != A2 & nchar(Minor) == 1 & !between(FRQ, 0.42, 0.58),
                 1 - FRQ, FRQ),
    OR = ifelse(Minor != A2 & nchar(Minor) == 1  & !between(FRQ, 0.42, 0.58),
                1 / OR, OR),
    A1 = ifelse(Minor != A2 & nchar(Minor) == 1  & !between(FRQ, 0.42, 0.58),
                A2, A1),
    A2 = ifelse(Minor != A2 & nchar(Minor) == 1  & !between(FRQ, 0.42, 0.58),
                Minor, A2)) %>%
  left_join(select(input$gencode, SNP, gencode_gene, gencode_dist, direction),
            by = "SNP") %>%
  left_join(select(input$gvc, SNP, gvc_gene, gvc_dist), by = "SNP") %>%
  left_join(input$snpeff, by = c("SNP" = "snp")) %>%
  mutate(
    ## replace Amino acid codes
    hgvs_p = str_replace(hgvs_p, "p.", ""),
    hgvs_p = str_replace_all(hgvs_p, "Ala", "A"),
    hgvs_p = str_replace_all(hgvs_p, "Cys", "C"),
    hgvs_p = str_replace_all(hgvs_p, "Asp", "D"),
    hgvs_p = str_replace_all(hgvs_p, "Glu", "E"),
    hgvs_p = str_replace_all(hgvs_p, "Phe", "F"),
    hgvs_p = str_replace_all(hgvs_p, "Gly", "G"),
    hgvs_p = str_replace_all(hgvs_p, "His", "H"),
    hgvs_p = str_replace_all(hgvs_p, "ILE", "I"),
    hgvs_p = str_replace_all(hgvs_p, "Lys", "K"),
    hgvs_p = str_replace_all(hgvs_p, "Leu", "L"),
    hgvs_p = str_replace_all(hgvs_p, "Met", "M"),
    hgvs_p = str_replace_all(hgvs_p, "Asn", "N"),
    hgvs_p = str_replace_all(hgvs_p, "Pro", "P"),
    hgvs_p = str_replace_all(hgvs_p, "Gln", "Q"),
    hgvs_p = str_replace_all(hgvs_p, "Arg", "R"),
    hgvs_p = str_replace_all(hgvs_p, "Ser", "S"),
    hgvs_p = str_replace_all(hgvs_p, "Thr", "T"),
    hgvs_p = str_replace_all(hgvs_p, "Val", "V"),
    hgvs_p = str_replace_all(hgvs_p, "Trp", "W"),
    hgvs_p = str_replace_all(hgvs_p, "Tyr", "Y"),
  ) %>%
  separate(hgvs_p, into = c("p1", "p2"),  sep = "\\d{1,4}", remove = FALSE) %>%
  mutate(
    ## update amino acids to use match minor allele
    aa = str_extract(hgvs_p,"\\d{1,4}"),
    hgvs_p_new = case_when(
      Minor == allele ~ glue("{p1}{aa}{p2}"),
      Minor != allele ~ glue("{p2}{aa}{p1}")
         )
  ) %>%
  select(-p1, -p2, -aa) %>%
  mutate(locus_ld = as.numeric(locus_ld))

# Export
message("\nExporting....", outfile, "\n")
write_csv(out, outfile)
