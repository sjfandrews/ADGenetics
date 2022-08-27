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
library(httr)

# Functions
## Calculated LD from LDlink
calc_ld <- function(snp_df) {
  message("Computing ld for chr ", snp_df %>% slice(1) %>% pull(chr), "...")

  if (nrow(snp_df) > 1) {
    out <- LDlinkR::LDmatrix(snps = snp_df$SNP,
                             pop = "EUR", r2d = "r2",
                             token = "dbae3c3dc2d4"
    )
  } else {
    out <- NULL
  }
  out
}

## Taging SNPs from LDlink
tag_snps <- function(snp_df) {
  message("Computing ld for chr ", snp_df %>% slice(1) %>% pull(chr), "...")

  if (nrow(snp_df) > 1) {
    out <- LDlinkR::SNPclip(snps = snp_df$SNP,
                            pop = "EUR",
                            r2_threshold = "0.1",
                            maf_threshold = "0.001",
                            token = "dbae3c3dc2d4"
    )
  } else {
    out <- NULL
  }
  out
}

# Import datasets
other <- read_csv("intermediate/other_loci.csv",
                  col_types = list(A2 = col_character()))
snps <- list(
  Lambert = "intermediate/lambert_loci.csv",
  Kunkle = "intermediate/kunkle_loci.csv",
  Marioni = "intermediate/marioni_loci.csv",
  Jansen = "intermediate/jansen_loci.csv",
  Wightman = "intermediate/wightman_loci.csv",
  Bellenguez = "intermediate/bellenguez_loci.csv") %>%
  map(read_csv) %>%
  map(select, SNP, CHR, BP, A1, A2, any_of(c("GENE")), FRQ, OR, P) %>%
  imap_dfr(~ mutate(.x, study = .y)) %>%
  bind_rows(other %>% filter(!str_detect(GENE, "APOE"))) %>%
  arrange(CHR, BP)

# Effective sample size
## https://github.com/neurogenomics/MungeSumstats/blob/master/R/compute_sample_size_neff.R
## meta data
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

write_csv(adgwas, "intermediate/ad_gwas_meta_data.csv")

## GWAS VCF
snp_list <- snps %>% 
  left_join(adgwas, by = "study") %>% 
  group_by(CHR, BP) %>%
  slice(which.max(neff)) %>%
  ungroup() %>%
  select(SNP, CHR, BP, A1, A2, OR, P, study) %>%
  rename(STUDY = study)

snp_list_munged <- MungeSumstats::format_sumstats(path = snp_list,
                                 ref_genome = "GRCh37",
                                 allele_flip_check = TRUE,
                                 allele_flip_drop = FALSE,
                                 bi_allelic_filter = FALSE,
                                 return_data = TRUE,
                                 log_folder = "data/MungeSumstats",
                                 force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

snp_list_out <- snp_list_munged %>%
  bind_rows(anti_join(snp_list, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP)

snp_list_out %>%
  MungeSumstats::write_sumstats("intermediate/adgwas_loci.vcf",
                                write_vcf = TRUE)

## Unique variants
snp_list <- distinct(snps, SNP)
write_delim(snp_list, "intermediate/adgwas_variant_list.txt",
            col_names = FALSE)

## Nearest Gene

### Gencode v40
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz
gencode_path <- "~/Downloads/gencode.v40lift37.annotation.gtf.gz"
gencode <- rtracklayer::import(gencode_path) %>%
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

loci_gene <- snps %>%
  rename(chr = CHR, pos = BP) %>%
  distinct(SNP, chr, pos) %>%
  left_join(
    purrr::map2_dfr(.$chr, .$pos, nearest_gene_df, gencode),
    by = c("chr", "pos")) %>% 
  rename(gencode_gene = gene_name, gencode_dist = dist)

## Extract MAF from NCBI DBSNP
dbsnp_info <- snp_list %>%
  # filter(!str_detect(SNP, "APOE")) %>%
  pull(SNP) %>%
  rsnps::ncbi_snp_query()

global_maf <- dbsnp_info %>%
  select(query, gene, maf_population) %>%
  unnest(cols = c(maf_population)) %>%
  filter(study == "dbGaP_PopFreq") %>%
  group_by(query) %>%
  slice(which.max(MAF)) %>%
  ungroup() %>%
  mutate(
    AF = MAF,
    MAF = ifelse(MAF > 0.5, 1 - MAF, MAF),
    Major = ifelse(AF > 0.5, Minor, ref_seq),
    Minor = ifelse(AF < 0.5, Minor, ref_seq),
  ) %>%
  relocate(Major, .before = Minor) %>%
  mutate(gene = na_if(gene, "")) %>%
  rename(dbsnp_gene = gene, SNP = query)

## MAF from gnomad v2.1.1
#system(paste(
#  "bcftools view -R input/adgwas_loci.vcf",
#  "-Ov -o output/test.vcf",
#  "reference/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz"))

gnomadaf_raw <- vcfR::read.vcfR("intermediate/gnomad_af.vcf")

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

## Locus Definiations
### Genomic region (Bellenguez et al 2022. loci definition)
file_loci_unmerged <- "intermediate/adgwas_loci.bed"
file_loci_merged <- "intermediate/adgwas_loci_merged.bed"

snp_list_out %>%
  mutate(
    chrom = paste0("chr", CHR),
    start = BP - 500000,
    end = BP + 500000,
    start = ifelse(start < 0, 1, start),
  ) %>%
  select(chrom, start, end, SNP, STUDY) %>%
  unite(snp, SNP, STUDY) %>%
  write_tsv(file_loci_unmerged, col_names = FALSE)

# system(glue("bedtools merge -i {file_loci_unmerged} > {file_loci_merged}"))

region_locus <- glue("bedtools merge -i {file_loci_unmerged}") %>%
  system(intern = TRUE) %>%
  I() %>%
  read_tsv(col_names = FALSE) %>%
  rename(chr = X1, start = X2, end = X3) %>%
  mutate(locus = glue("{chr}:{start}-{end}"),
         chr = str_replace(chr, "chr", ""),
         chr = as.numeric(chr))

### LD Structure
ld <- snps %>%
  select(SNP, CHR, BP, A1, A2, study) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  mutate(chr = CHR) %>%
  nest(data = c(SNP, chr, BP, A1, A2, study)) %>%
  arrange(CHR) %>%
  mutate(ld = map(data, calc_ld),
         tags = map(data, tag_snps))

tags_prune <- ld %>%
  select(tags) %>%
  unnest(cols = tags) %>%
  mutate(tag = str_extract(Details, "rs[:digit:]+"),
         tag = ifelse(is.na(tag), RS_Number, tag)) %>%
  separate(Position, c("CHR", "POS")) %>%
  full_join(ld %>% select(data) %>% unnest(data) %>% select(SNP, chr, BP),
            by = c("RS_Number" = "SNP")) %>%
  mutate(CHR = str_replace(CHR, "chr", ""),
         CHR = as.numeric(CHR),
         POS = as.numeric(POS),
         # Fill in infromation from variants not in 1KG
         CHR = ifelse(is.na(CHR), chr, CHR),
         POS = ifelse(is.na(POS), BP, POS),
         tag = ifelse(is.na(tag), RS_Number, tag),
         tag = as_factor(tag)) %>%
  select(-chr, -BP) %>%
  arrange(CHR, POS) %>%
  mutate(
    locus_ld = lvls_revalue(tag, as.character(seq_along(fct_unique(tag))))
  )

tags <- snps %>%
  # get the minimum P value for each variant
  arrange(CHR, BP, P) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  select(SNP, P) %>%
  # Join the minimum P value for each variant onto the tag df
  left_join(tags_prune, ., by = c("RS_Number" = "SNP")) %>%
  # Find the minimum P value for each ld locus
  group_by(locus_ld) %>%
  arrange(P) %>%
  distinct(locus_ld, .keep_all = TRUE) %>%
  ungroup %>%
  # Get the tagging SNP by clump and join back onto the tag df
  select(tag_clump = RS_Number, locus_ld) %>%
  left_join(tags_prune, ., by = "locus_ld")

### Annotate with cytogenic bands
# Rtracklayer downloads different cgbands than a command line query
# bands <-
#  rtracklayer::browserSession() %>%
#  rtracklayer::ucscTableQuery(table = "cytoBand", genome = "hg19") %>%
#  rtracklayer::getTable() %>%
#  tibble::as_tibble()


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

### Merge regions, LD, and cytogenic bands

loci <- tags %>%
  mutate(snp_start = POS, snp_end = POS + 1) %>%
  fuzzyjoin::genome_left_join(region_locus,
    by = c("CHR" = "chr", "snp_start" = "start", "snp_end" = "end")) %>%
  fuzzyjoin::genome_left_join(bands,
    by = c("CHR" = "chr", "snp_start" = "chromStart",
           "snp_end" = "chromEnd")) %>%
  select(SNP = RS_Number, POS, Alleles, Details,
         tag, tag_clump, locus_ld, locus, cytoband)

## snpEff
file_ann <- "intermediate/adgwas_loci.ann.vcf"

"java -Xmx8g -jar snpEff.jar GRCh37.75 %s > %s" %>%
  sprintf("intermediate/adgwas_loci.vcf", file_ann) %>%
  system()

raw_snpeff <- file_ann %>%
  vcfR::read.vcfR() %>%
  vcfR::vcfR2tidy(., info_only = FALSE)

ann_cols <- raw_snpeff$meta %>%
  filter(Tag == "INFO" & ID == "ANN") %>%
  pull(Description) %>%
  str_match("'(.+)'") %>%
  .[[2]] %>%
  str_split(" \\| ") %>%
  .[[1]]

ann <- raw_snpeff$gt %>%
  select(ChromKey, POS, SNP = gt_SNP) %>%
  left_join(., raw_snpeff$fix, by = c("ChromKey", "POS")) %>%
  select(-ChromKey, -ID, -QUAL, -FILTER, -LOF, -NMD) %>%
  mutate(ANN = strsplit(ANN, ",")) %>%
  unnest(cols = ANN) %>%
  separate(ANN, into = ann_cols, sep = "\\|") %>%
  janitor::clean_names() %>%
  mutate(chrom = as.numeric(chrom))

snpeff <- ann %>%
  filter(transcript_bio_type %in% c("", "protein_coding")) %>%
  group_by(snp) %>%
  slice_head() %>%
  ungroup() %>%
  arrange(chrom, pos)

## AD Gene verification committee
# https://adsp.niagads.org/index.php/gvc-top-hits-list/
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

gvc_genes_raw <- read_sheet(gvc_url, sheet = 2, skip = 1)

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


## Join datasets together

out <- other %>%
  filter(str_detect(GENE, "APOE")) %>%
  bind_rows(snps, .) %>%
  left_join(select(global_maf, SNP, dbsnp_gene, Major, Minor, MAF, AF),
            by = "SNP") %>%
  left_join(select(gnomad_maf, SNP, gnomad_minor, gnomad_af, gnomad_maf),
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
  left_join(select(loci_gene, SNP, gencode_gene, gencode_dist, direction),
            by = "SNP") %>%
  left_join(select(loci_gvc, SNP, gvc_gene, gvc_dist), by = "SNP") %>%
  left_join(snpeff, by = c("SNP" = "snp")) %>%
  mutate(locus_ld = as.numeric(locus_ld)) #%>%
  #mutate(locus_ld = ifelse(locus_ld >= 88, locus_ld + 1, locus_ld))

write_csv(out, "intermediate/ad_loci.csv")
write_rds(ld, "intermediate/ad_ld.rds")