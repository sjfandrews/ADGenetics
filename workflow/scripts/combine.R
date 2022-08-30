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

adgwas.path = snakemake@input[['adgwas']]
other.path = snakemake@input[['other']]
tags.path = snakemake@input[['tags']]
locus.path = snakemake@input[['locus']]
bands.path = snakemake@input[['bands']]
dbsnp.path =  snakemake@input[['dbsnp']]
gnomad.path = snakemake@input[['gnomad']]
gencode.path = snakemake@input[['gencode']]
gvc.path = snakemake@input[['gvc']]
snpeff.path = snakemake@input[['snpeff']]

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

### For updating APOE locus 
apoe_locus <- loci %>% filter(SNP == "rs429358") %>% slice(1) %>% select(locus, cytoband, locus_ld)


## Join datasets together
out <- other %>% 
  filter(str_detect(GENE, "APOE")) %>%
  bind_rows(snps, .) %>% 
  left_join(select(dbsnp, SNP, dbsnp_gene, Major, Minor, MAF, AF),
            by = "SNP") %>%
  left_join(select(gnomad, SNP, gnomad_minor, gnomad_af, gnomad_maf),
            by = "SNP") %>%
  dplyr::rename(global_maf = MAF, global_af = AF) %>% 
  left_join(select(loci, SNP, locus, cytoband, locus_ld), by = c("SNP")) %>%
  relocate(locus, cytoband, locus_ld) %>%
  mutate(
    ## Manual fixing missing data for rs139643391 & rs149080927
    A1 = ifelse(SNP == "rs139643391", "TC", A1),
    A2 = ifelse(SNP == "rs139643391", "T", A2),
    FRQ = ifelse(SNP == "rs139643391", 0.131, FRQ),
    OR = ifelse(SNP == "rs139643391", 0.94, OR),
    global_maf = ifelse(SNP == "rs139643391", 0.087446, global_maf),
    gnomad_maf = ifelse(SNP == "rs149080927", 0.3829, gnomad_maf),
  ) %>%
  mutate(
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
  left_join(select(gencode, SNP, gencode_gene, gencode_dist, direction),
            by = "SNP") %>%
  left_join(select(gvc, SNP, gvc_gene, gvc_dist), by = "SNP") %>%
  left_join(snpeff, by = c("SNP" = "snp")) %>%
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
  separate(hgvs_p, into = c("p1", "p2"),  sep = "\\d{1,4}", remove = F) %>%
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































