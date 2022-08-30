library(MungeSumstats)
library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(plyr)
library(dplyr)
library(stringr)
library(forcats)
library(tabulizer)
library(fuzzyjoin)

data(sumstatsColHeaders)

# Snakemake input & output
lambert.path <- snakemake@input[['lambert']]
kunkle.path <- snakemake@input[['kunkle']]
marioni.path <- snakemake@input[['marioni']]
jansen.path <- snakemake@input[['jansen_table']]
jansen_ss.path <- snakemake@input[['jansen_ss']]
wightman.path <- snakemake@input[['wightman_table']]
wightman_ss.path <- snakemake@input[['wightman_ss']]
bellenguez.path <- snakemake@input[['bellenguez']]

lambert_out.path = snakemake@output[['lambert_out']]
kunkle_out.path = snakemake@output[['kunkle_out']]
marioni_out.path = snakemake@output[['marioni_out']]
jansen_out.path = snakemake@output[['jansen_out']]
wightman_out.path = snakemake@output[['wightman_out']]
bellenguez_out.path = snakemake@output[['bellenguez_out']]
other_out.path = snakemake@output[['other_out']]

## Misc functions
calc_se_from_beta_and_l95 = function(or, l95){
  logOR = log(or)
  logOR_L95 = log(l95)
  logOR_SE = (logOR-logOR_L95)/1.96
  logOR_SE
}

calc_beta_from_z_n_freq = function(Z, FRQ, N){
  out = Z/sqrt(2*FRQ*(1-FRQ)*(N+Z^2))
  out
}

## Lambert 2013
message("\nMunging Lambert 2013\n")

lambert_raw <- extract_tables(lambert.path)
lambert_clean <- as_tibble(lambert_raw[[1]]) %>% 
  slice(-c(1,2,3,4,23,34)) %>%
  select(V1:V5, V12:V13) %>%
  magrittr::set_colnames(c("snp", "chr", "pos", "gene", "alleles", "or", "p", "hetpval")) %>%
  separate(alleles, c("alleles", "maf"), sep = " ") %>%
  mutate(snp = na_if(snp, ""), 
         effect = rep(c("or", "ci"), 21)) %>%
  fill(snp, .direction = "down") %>%
  pivot_wider(names_from = effect, values_from = or) %>%
  fill(ci, .direction = "up") %>%
  filter(!is.na(or)) %>%
  mutate(ci = str_replace_all(ci, "\\(|\\)", "")) %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("otherAllele", "effectAllele"), sep = "/") %>%
  mutate(p = str_replace_all(p, " × 10−", "e-"), 
         maf = as.numeric(maf),
         p = as.numeric(p),
         or = as.numeric(or), 
         l95 = as.numeric(l95), 
         u95 = as.numeric(u95), 
         beta = log(or), 
         se = calc_se_from_beta_and_l95(or, l95),
         chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         snp = str_replace(snp, "rs3865444g", "rs3865444"), 
         snp = str_replace(snp, "rs8093731g", "rs8093731")) %>%
  read_sumstats( .,
    standardise_headers = T,
    mapping_file = sumstatsColHeaders
  ) %>%
  filter(P < 5e-8) %>%
  mutate(CHR = as.numeric(CHR)) %>%
  arrange(CHR, BP) %>%
  as_tibble()

lambert <- MungeSumstats::format_sumstats(path=lambert_clean, 
                               ref_genome="GRCh37", 
                               allele_flip_check = FALSE, 
                               allele_flip_drop = FALSE, 
                               bi_allelic_filter = FALSE,
                               return_data = TRUE, 
                               log_folder = "data/MungeSumstats", 
                               force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(lambert_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

## kunkle 2019
message("\nMunging Kunkle 2019\n")

kunkle_raw <- extract_tables(kunkle.path)
k1 <- as_tibble(kunkle_raw[[1]]) %>% 
  slice(-c(1:7,27,30)) %>%
  select(-V3) %>%
  mutate(V1 = str_replace_all(V1, "!×!10−| ×!10−|x 10−", "e-"), 
         V2 = str_replace_all(V2, "!×!10−| ×!10−|x 10−", "e-")) %>%
  separate(V1, c("snp", 'chr', 'pos', 'gene','alleles', 'maf', 'or1', 'ci1', 'p1', 
                 'or2', 'ci2', 'p2'), sep = " ") %>%
  separate(V2, c('or', 'ci', 'p'), sep = " ") %>% 
  mutate(or = na_if(or, ""),
         or = ifelse(is.na(or), or1, or), 
         ci = ifelse(is.na(ci), ci1, ci), 
         p = ifelse(is.na(p), p1, p)) %>%
  select(snp, chr, pos, gene, alleles, maf, or, ci, p) %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("otherAllele", "effectAllele"), sep = "/")

k2 <- as_tibble(kunkle_raw[[2]]) %>% 
  slice(-c(1:4,14:16)) %>% 
  select(V1, V2, V3, V4, V6, V11, V12) %>%
  mutate(V12 = str_replace_all(V12, "!×!10−| ×!10−|x 10−", "e-")) %>%
  magrittr::set_colnames(c("snp", "chr", "pos", "gene", "maf", "or", "p")) %>%
  separate(gene, c("gene",'alleles'), sep = " ") %>%
  separate(or, c("or",'ci'), sep = " ") %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("otherAllele", "effectAllele"), sep = "/")

kunkle_clean <- bind_rows(
  k1, k2
) %>%
  mutate(p = as.numeric(p),
         or = as.numeric(or), 
         l95 = as.numeric(l95), 
         u95 = as.numeric(u95), 
         beta = log(or), 
         maf = as.numeric(maf), 
         se = calc_se_from_beta_and_l95(or, l95),
         chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         snp = str_replace(snp, "rs7920721g", "rs7920721")) %>%
  arrange(chr, pos) %>% 
  filter(p < 5e-8) %>%
  group_by(snp) %>%
  slice(which.min(p)) %>% 
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

kunkle <- MungeSumstats::format_sumstats(path=kunkle_clean, 
                                          ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                          log_folder = "data/MungeSumstats", 
                                          force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(kunkle_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

## Marioni 2018
message("\nMunging Marioni 2018\n")

marioni_ss <- read_tsv(marioni.path) %>%
  select(SNP = DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, OR, SE, Z, P) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

marioni <- MungeSumstats::format_sumstats(path=marioni_ss, 
                                         ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(marioni_ss, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) 

## Jansen 2019 
message("\nMunging Jansen 2019\n")

jansen_raw <- extract_tables(jansen.path)
jansen_ss <- read_tsv(jansen_ss.path) %>%
  select(-uniqID.a1a2, Nsum, -Neff, -dir, -Z) %>%
  rename(effectAllele = A1, otherAllele = A2)

jansen_clean <- as_tibble(jansen_raw[[1]]) %>%
  slice(-c(1:5, 19)) %>%
  mutate(V1 = str_replace_all(V1, " 2.59 ×", " 2.59E10-4"),
         V1 = str_replace_all(V1, "!×!| ×!|x |!× !|!× | × ", "E"), 
         V1 = str_replace_all(V1, " NA ", " NA NA "), 
         V1 = str_replace_all(V1, "− ", "-")) %>%
  separate(V1, into = c("locus", "chr", "gene", "snp1", "P1", "SNP2", "P2", "snp", "pos", "a2", "a1", "maf", "z", "p", "dir"), sep = " ") %>%
  select(snp, gene) %>%
  left_join(jansen_ss, by = c("snp" = "SNP")) %>%
  rename(n = Nsum) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  filter(P < 5e-8)

jansen <- MungeSumstats::format_sumstats(path=jansen_clean, 
                                         ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(jansen_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE, OR = exp(BETA))

print(jansen, n = Inf)

## Wightman 2021
message("\nMunging Wightman 2021\n")

wightman_raw <- extract_tables(wightman.path)
wightman_ss <- read_tsv(wightman_ss.path) %>%
  select(-p, -N) %>%
  rename(effectAllele = testedAllele)

wightman_clean <- as_tibble(wightman_raw[[1]]) %>%
  slice(-c(1:3)) %>%
  mutate(V1 = str_replace_all(V1, "!×!10−| ×!10−|x 10−|!× !10−|!× 10−| × 10−", "e-"), 
         V1 = str_replace_all(V1, ",", ""), 
         V2 = str_replace_all(V2, ",", "")) %>%
  separate(V1, into = c("locus", "gene", "pos", "snp", "ea", "freq", "pval"), sep = " ") %>%
  separate(pos, into = c("chr", "pos"), sep = ":") %>%
  rename(n = V2) %>%
  # select(snp, chr, pos, gene, freq, p, n)  %>%
  mutate(chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         freq = as.numeric(freq),
         n = as.numeric(n), 
         pval = as.numeric(pval), 
         pval = ifelse(is.na(pval), 1.0e-300, pval)) %>%
  left_join(wightman_ss, by = c("chr", "pos" = "PosGRCh37")) %>%
  mutate(freq = ifelse(ea != effectAllele, 1-freq, freq)) %>% 
  select(-ea) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) 

wightman <- MungeSumstats::format_sumstats(path=wightman_clean, 
                                         ref_genome="GRCh37", 
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(wightman_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = ifelse(is.infinite(Z), 500, Z),
         BETA = calc_beta_from_z_n_freq(Z, FRQ, N), 
         OR = exp(BETA), 
         )

print(wightman, n = Inf)

## Bellenguez 2022 
message("\nMunging Bellenguez 2022\n")
### Adding APOE locus added manually from Kunkle summary statistics 

bellenguez_raw <- extract_tables(bellenguez.path)

bellenguez_clean <- bind_rows(
  bellenguez_raw[[1]] %>% 
    as_tibble() %>%
    slice(4:nrow(.)) %>%
    separate(col = "V1", into = c("snp", "chr", "pos", "gene", 'locus', 'alleles', 'maf', 'or', 'ci', 'p'), 
             sep = " "), 
  
  bellenguez_raw[[2]] %>% 
    as_tibble() %>%
    slice(6:nrow(.)) %>%
    mutate(V1 = str_replace(V1, "IGH gene cluster", "IGH_gene_cluster")) %>%
    separate(col = "V1", into = c('locus', "snp", "chr", "pos", "gene", 'alleles', 'maf', 'or', 'ci', 'p'), 
             sep = " ")
) %>%
  separate(ci, c("l95", "u95"), sep = "–") %>%
  separate(alleles, c("effectAllele", "otherAllele"), sep = "/") %>%
  mutate(p = str_replace_all(p, "!×!10−| ×!10−|x 10−|!× !10−|!× 10−| × 10−", "e-"),
         or = as.numeric(or), 
         l95 = as.numeric(l95), 
         u95 = as.numeric(u95), 
         p = as.numeric(p), 
         chr = as.numeric(chr), 
         pos = as.numeric(pos), 
         beta = log(or), 
         se = calc_se_from_beta_and_l95(or, l95)) %>%
  add_row(snp = "rs429358", chr = 19, pos = 44908684, gene = "APOE", locus = "APOE", 
          effectAllele = "C", otherAllele = "T", maf = '0.216', 
          or = 3.32, l95 = 3.20, u95 = 3.45, p = 1.2e-881, beta = 1.199, se = 0.019) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  arrange(CHR, BP)
  
bellenguez_munged <- MungeSumstats::format_sumstats(path=bellenguez_clean, 
                                         ref_genome="GRCh38", 
                                         convert_ref_genome = "GRCh37",
                                         allele_flip_check = FALSE, 
                                         allele_flip_drop = FALSE, 
                                         bi_allelic_filter = FALSE,
                                         return_data = TRUE, 
                                         log_folder = "data/MungeSumstats", 
                                         force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR))

### update BP for SNPs not lifted over  
bellenguez <- bellenguez_munged %>% 
  bind_rows(anti_join(bellenguez_clean, ., by = "SNP")) %>%
  mutate(BP = case_when(
           SNP == "rs139643391" ~ 203743439,
           SNP == "rs6605556" ~ 32583099,
           SNP == "rs1160871" ~ 28168746,
           SNP == "rs35048651" ~ 1631341,
           SNP == "rs199515" ~ 44856641,
           SNP == "rs149080927" ~ 1854258,
           SNP == "rs587709" ~ 54771451,
           TRUE ~ BP
         )) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

print(bellenguez, n = Inf)

## Adding other variants of interst
other <- tribble(
  ~SNP, ~CHR, ~BP, ~A1, ~A2, ~GENE, ~FRQ, ~OR, ~P, ~study,
  "rs63750847", 21, 27269932, "C", "T", "APP", 0.0003, 0.23, 4.19e-5, "Jonsson",
  "APOE 22", 19, 45412079, "TT/CC", "TT/TT", "APOE ε2/ε2", 0.007, 0.42, 7.1e-4, "Reiman",
  "APOE 23", 19, 45412079, "TT/CC", "TT/CT", "APOE ε2/ε3", 0.126, 0.58, 3.1e-25, "Reiman",
  "APOE 24", 19, 45412079, "TT/CC", "CT/CT", "APOE ε2/ε4", 0.027, 2.49, 9.8e-32, "Reiman",
  "APOE 34", 19, 45412079, "TT/CC", "CT/CC", "APOE ε3/ε4", 0.247, 3.78, 1e-300, "Reiman",
  "APOE 44", 19, 45412079, "TT/CC", "CC/CC", "APOE ε4/ε4", 0.028, 12.28, 1e-300, "Reiman",
)

## Export Cleaned datasets 
write_csv(lambert, lambert_out.path)
write_csv(kunkle, kunkle_out.path)
write_csv(marioni, marioni_out.path)
write_csv(jansen, jansen_out.path)
write_csv(wightman, wightman_out.path)
write_csv(bellenguez, bellenguez_out.path)
write_csv(other, other_out.path)
































































