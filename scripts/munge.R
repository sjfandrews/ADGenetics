library(MungeSumstats)
library(tidyverse)
library(tabulizer)
library(fuzzyjoin)

data(sumstatsColHeaders)
setwd("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/")

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
lambert_raw <- extract_tables("resources/Lambert2013/www.nature.com 11_22_2019, 11_45_10 AM.pdf")
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
  filter(p < 5e-8) %>%
  arrange(chr, pos) %>%
  read_sumstats( .,
    standardise_headers = T,
    mapping_file = sumstatsColHeaders
  ) %>%
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
kunkle_raw <- extract_tables("resources/www.nature.com 11_18_2019, 11_19_16 AM.pdf")

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
  as_tibble() 

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

print(kunkle, n = Inf)

## Marioni 2018
marioni_ss <- read_tsv("resources/Marioni2018_gws.txt") %>%
  select(SNP = DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, OR, SE, Z, P) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble()

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
jansen_raw <- extract_tables("resources/www.nature.com 11_25_2019, 7_56_44 AM.pdf")
jansen_ss <- read_tsv("resources/Jansen2019_gws.txt") %>%
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
wightman_raw <- extract_tables("resources/www.nature.com 02_10_2021, 15_58_50.pdf")
wightman_ss <- read_tsv("resources/wightman_gws.txt") %>%
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
  as_tibble() 

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
bellenguez_raw <- extract_tables("~/Dropbox/Research/PostDoc-MSSM/Neurogenomics/resources/Bellenguez2022.pdf")

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
  arrange(chr, pos) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble()
  
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


## Naj et al 2021
naj_raw <- extract_tables("resources/Naj_et_al-2021-medRxiv.pdf")

naj_t2 <- naj_raw[[2]] %>%
  as_tibble() %>%
  slice(4:nrow(.)) %>%
  magrittr::set_colnames(c("snp", "chr", "pos", "locus", "effectAllele", "otherAllele", 
                         "stage2",  
                         "ukb_eaf", "ukb_or", "ukb_p", 
                         "stage3")) %>%
  select(-stage2, -ukb_eaf, -ukb_or, -ukb_p) %>%
  separate(stage3, c("eaf", "or", 'l95', 'u95', 'p'), sep = " ") %>%
  mutate(
    l95 = as.numeric(str_replace_all(l95, "\\(|,", "")), 
    u95 = as.numeric(str_replace_all(u95, "\\)", "")), 
    or = as.numeric(or), 
    eaf = as.numeric(eaf), 
    chr = as.numeric(chr), 
    pos = as.numeric(pos), 
    p = as.numeric(str_replace_all(p, "E-", "e-")), 
    )

naj_t4 <- naj_raw[[4]] %>%
  as_tibble() %>%
  slice(4:nrow(.)) %>%
  magrittr::set_colnames(c("snp", "chr", "pos", "locus", "effectAllele", "otherAllele", 
                           "stage2",  
                           "ukb_eaf", "ukb_or", "ukb_p", 
                           "stage3")) %>%
  select(-stage2, -ukb_eaf, -ukb_or, -ukb_p) %>%
  separate(stage3, c("eaf", "or", 'l95', 'u95', 'p'), sep = " ") %>%
  mutate(
    l95 = as.numeric(str_replace_all(l95, "\\(|,", "")), 
    u95 = as.numeric(str_replace_all(u95, "\\)", "")), 
    or = as.numeric(or), 
    eaf = as.numeric(eaf), 
    chr = as.numeric(chr), 
    pos = as.numeric(pos), 
    p = as.numeric(str_replace_all(p, "E-", "e-")), 
  )

naj_clean <- bind_rows(
  naj_t2, naj_t4
) %>%
  arrange(chr, pos) %>%
  filter(p < 5e-8) %>%
  mutate(beta = log(or), 
         se = calc_se_from_beta_and_l95(or, l95)) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() 

naj <- MungeSumstats::format_sumstats(path=naj_clean, 
                                             ref_genome="GRCh37", 
                                             allele_flip_check = FALSE, 
                                             allele_flip_drop = FALSE, 
                                             bi_allelic_filter = FALSE,
                                             return_data = TRUE, 
                                             log_folder = "data/MungeSumstats", 
                                             force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(naj_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(Z = BETA / SE)

## Bis 2017 

bis_raw <- extract_tables("resources/Bis2017.pdf")

bis_raw[[1]] %>%
  as_tibble() %>%
  slice(4:nrow(.)) %>%
  magrittr::set_colnames(c("snp", "gene", "discovery", "modles", "replication", 'all1', "all2")) %>%
  separate(snp, c("variant", "snp"), sep = " ") %>%
  separate(variant, c("chr", "pos", "a1", "a2"), sep = ":")

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
write_csv(lambert, "intermediate/lambert_loci.csv")
write_csv(kunkle, "intermediate/kunkle_loci.csv")
write_csv(marioni, "intermediate/marioni_loci.csv")
write_csv(jansen, "intermediate/jansen_loci.csv")
write_csv(wightman, "intermediate/wightman_loci.csv")
write_csv(bellenguez, "intermediate/bellenguez_loci.csv")
write_csv(other, "intermediate/other_loci.csv")


































































