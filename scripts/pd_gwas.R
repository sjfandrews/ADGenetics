library(MungeSumstats)
library(tidyverse)
library(tabulizer)
library(fuzzyjoin)
library(plotgardener)

data(sumstatsColHeaders)
setwd("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/")

#Chang 2017
# chang_raw <- extract_tables("resources/Chang2017.pdf") 
# 
# chang_munged <- bind_rows(
#   chang_raw %>%
#     magrittr::extract2(1) %>%
#     as_tibble() %>%
#     slice(-1) %>%
#     select(-V6, -V7, -V8, -V9, -V10, -V11, -V12, -V13, -V14, -V15, -V16) %>%
#     magrittr::set_colnames(c('pos', 'rsid', 'gene', 'allele', 'eaf', 'p', 'or', 'ci')) %>%
#     mutate(gene = case_when(
#       rsid == "rs199347" ~ "KLHL7, NUPL2, GPNMB", 
#       rsid == "rs17649553" ~ "ARHGAP27, CRHR1, SPPL2C, MAPT, STH, KANSL1", 
#       TRUE ~ gene
#     )) %>%
#     filter(rsid != "") %>%
#     separate(allele, c('effect_allele', "other_allele"), sep = "/") %>% 
#     separate(ci, c('lci', "uci"), sep = "–") %>% 
#     separate(pos, c('chr', "pos"), sep = ":") %>% 
#     mutate(
#       rsid = str_extract(rsid, "rs[:digit:]*"),
#       chr = as.numeric(chr),
#       pos = str_extract(pos, "[:digit:]*"),
#       pos = as.numeric(pos),
#       eaf = as.numeric(eaf),
#       or = as.numeric(or),
#       lci = as.numeric(lci),
#       uci = as.numeric(uci),
#       p = str_replace_all(p, " × 10−", "e-"),
#       p = as.numeric(p)
#     ) , 
#   
#   chang_raw %>%
#     magrittr::extract2(2) %>%
#     as_tibble() %>%
#     slice(4:n()) %>%
#     select(-V6, -V7, -V8, -V9) %>%
#     magrittr::set_colnames(c('pos', 'rsid', 'gene', 'allele', 'eaf', 'p', 'or', 'ci')) %>%
#     mutate(
#       gene = case_when(
#         rsid == "rs12497850" ~ "NCKIPSD, CDC71",
#         rsid == "rs143918452" ~ "ALAS1, TLR9, DNAH1, BAP1, PHF7, NISCH, STAB1, ITIH3, ITIH4",
#         rsid == "rs78738012" ~ "ANK2, CAMK2D",
#         rsid == "rs2280104" ~ "SORBS3, PDLIM2, C8orf58, BIN3",
#         rsid == "rs601999" ~ "ATP6V0A1, PSMC3IP, TUBG2",
#         TRUE ~ gene
#       )
#     ) %>%
#     filter(rsid != "") %>%
#     separate(allele, c('effect_allele', "other_allele"), sep = "/") %>% 
#     separate(ci, c('lci', "uci"), sep = "–") %>% 
#     separate(pos, c('chr', "pos"), sep = ":") %>% 
#     mutate(
#       rsid = str_extract(rsid, "rs[:digit:]*"),
#       chr = as.numeric(chr),
#       pos = str_extract(pos, "[:digit:]*"),
#       pos = as.numeric(pos),
#       eaf = str_extract(eaf, "0.[:digit:]*"),
#       eaf = as.numeric(eaf),
#       or = as.numeric(or),
#       lci = as.numeric(lci),
#       uci = as.numeric(uci),
#       p = str_replace_all(p, " × 10−", "e-"),
#       p = as.numeric(p), 
#       p = ifelse(rsid == "rs601999", 8.03e-9 , p),
#       or = ifelse(rsid == "rs601999", 0.93 , or),
#       lci = ifelse(rsid == "rs601999", 0.915 , lci),
#       uci = ifelse(rsid == "rs601999", 0.941 , uci),
#     ) 
# ) %>%
#   filter(p < 5e-8) %>%
#   arrange(chr, pos) %>%
#   filter(!is.na(lci)) %>%
#   read_sumstats( .,
#                  standardise_headers = T,
#                  mapping_file = sumstatsColHeaders
#   ) %>%
#   as_tibble()
# 
# # Nalls 2019
# ## Table 1
# nalls_raw <- extract_tables("resources/Nalls2019.pdf") 
# nalls_munged <- nalls_raw %>%
#   magrittr::extract2(1) %>%
#   as_tibble() %>%
#   slice(4:n()) %>%
#   select(-V4, -V7, -V9, -V11, -V14, -V16, -V18, -V19, -V20, -V21, -V22, -V23, -V24, -V25) %>%
#   magrittr::set_colnames(c('rsid', 'chr', 'bp', 'gene', 'effect_allele', 'other_allele', 'eaf', 'or', 'b', 'se', 'p')) %>%
#   separate(or, c('or', "ci"), sep = " ") %>% 
#   separate(ci, c('lci', "uci"), sep = "–") %>%
#   mutate(
#     chr = as.numeric(chr),
#     eaf = str_replace_all(eaf, "·", "."), 
#     eaf = as.numeric(eaf),
#     or = str_replace_all(or, "·", "."), 
#     or = as.numeric(or),
#     lci = str_replace_all(lci, "\\(", ""), 
#     lci = str_replace_all(lci, "·", "."), 
#     lci = as.numeric(lci),
#     uci = str_replace_all(uci, "\\)", ""),
#     uci = str_replace_all(uci, "·", "."), 
#     uci = as.numeric(uci),
#     b = str_replace_all(b, "–", "-"), 
#     b = str_replace_all(b, "·", "."), 
#     b = as.numeric(b),
#     se = str_replace_all(se, "·", "."), 
#     se = as.numeric(se),
#     p = str_replace_all(p, " × 10–", "e-"),
#     p = str_replace_all(p, "·", "."), 
#     p = as.numeric(p)
#   ) %>%
#   filter(p < 5e-8) %>%
#   arrange(chr, bp) %>%
#   read_sumstats( .,
#                  standardise_headers = T,
#                  mapping_file = sumstatsColHeaders
#   ) %>%
#   as_tibble() %>%
#   mutate(CHR = as.numeric(CHR))
# 
# 
# nalls <- MungeSumstats::format_sumstats(path=nalls_munged, 
#                                         ref_genome="GRCh37", 
#                                         allele_flip_check = FALSE, 
#                                         allele_flip_drop = FALSE, 
#                                         bi_allelic_filter = FALSE,
#                                         return_data = TRUE, 
#                                         log_folder = "data/MungeSumstats", 
#                                         force_new = TRUE) %>%
#   as_tibble() %>%
#   mutate(CHR = as.numeric(CHR)) %>%
#   bind_rows(anti_join(nalls_munged, ., by = c("CHR", "BP"))) %>%
#   arrange(CHR, BP) %>%
#   mutate(Z = BETA / SE)

## Supplementary Table 2
nalls_raw <- readxl::read_xlsx("resources/Nalls_ST2.xlsx") %>% 
  janitor::clean_names()
nalls_smr_raw <- readxl::read_xlsx("resources/Nalls_ST6.xlsx") %>% 
  janitor::clean_names()

nalls_clean <- nalls_raw %>%
  filter(failed_final_filtering_and_qc == 0) %>%
  select(locus_number, snp, chr, bp, nearest_gene, qtl_nominated_gene_nearest_qtl, 
         effect_allele, other_allele, effect_allele_frequency, beta_all_studies, 
         se_all_studies, p_all_studies) %>%
  rename(locus = locus_number, b = beta_all_studies, se = se_all_studies, p = p_all_studies) %>%
  read_sumstats( .,
                 standardise_headers = T,
                 mapping_file = sumstatsColHeaders
  ) %>%
  as_tibble() %>%
  mutate(
    LOCUS = as.numeric(LOCUS),
    CHR = as.numeric(CHR),
    A1 = toupper(A1), 
    A2 = toupper(A2), 
    QTL_NOMINATED_GENE_NEAREST_QTL = ifelse(is.na(QTL_NOMINATED_GENE_NEAREST_QTL), "null", QTL_NOMINATED_GENE_NEAREST_QTL)
  ) %>%
  arrange(CHR, BP)


nalls_munged <- MungeSumstats::format_sumstats(path=nalls_clean, 
                                        ref_genome="GRCh37", 
                                        allele_flip_check = TRUE, 
                                        allele_flip_drop = FALSE, 
                                        bi_allelic_filter = TRUE,
                                        return_data = TRUE, 
                                        log_folder = "data/MungeSumstats", 
                                        force_new = TRUE) %>%
  as_tibble() %>%
  mutate(CHR = as.numeric(CHR)) %>%
  bind_rows(anti_join(nalls_clean, ., by = c("CHR", "BP"))) %>%
  arrange(CHR, BP) %>%
  mutate(
    A3 = A1,
    A1 = case_when(
      SNP %in% c("rs76763715", "rs6825004", "rs1867598", "rs8087969") ~ A2,
      TRUE ~ A1
    ), 
    A2 = case_when(
      SNP %in% c("rs76763715", "rs6825004", "rs1867598", "rs8087969") ~ A3,
      TRUE ~ A2
    ), 
    FRQ = case_when(
      SNP %in% c("rs76763715", "rs6825004", "rs1867598", "rs8087969") ~ 1-FRQ,
      TRUE ~ FRQ
    ), 
    BETA = case_when(
      SNP %in% c("rs76763715", "rs6825004", "rs1867598", "rs8087969") ~ BETA * -1,
      TRUE ~ BETA
    )
  ) %>% 
  select(-A3) %>%
  mutate(Z = BETA / SE, 
         QTL_NOMINATED_GENE_NEAREST_QTL = ifelse(QTL_NOMINATED_GENE_NEAREST_QTL == "null", NA, QTL_NOMINATED_GENE_NEAREST_QTL), 
         MAF = ifelse(FRQ > 0.5, 1 - FRQ, FRQ), 
         BETA = ifelse(FRQ > 0.5, BETA * -1, BETA), 
         OR = exp(BETA)
         ) 



## Export
## GWAS VCF
nalls_munged %>%
  select(SNP, CHR, BP, A1, A2, FRQ, BETA, SE, P, OR) %>%
  data.table::as.data.table() %>%
  MungeSumstats::write_sumstats(., "intermediate/pdgwas_loci.vcf", write_vcf =T)

### MAF from gnomad v2.1.1
system('bcftools view -R input/pdgwas_loci.vcf -Ov -o output/pdgwas_loci_maf.vcf reference/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz')

pd_gnomadaf.raw <- vcfR::read.vcfR('~/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/pdgwas_loci_maf.vcf') 

pd_gnomad_af <- pd_gnomadaf.raw %>%
  vcfR::vcfR2tidy(., info_only = T) %>%
  magrittr::extract2("fix") %>%
  select(CHROM, POS, ID, REF, ALT, AF, AF_raw, AF_afr, AF_eas, AF_amr, AF_nfe, AF_asj, AF_fin) %>%
  mutate(AF = as.numeric(AF), 
         gnomad_maf = ifelse(AF > 0.5, 1 - AF, AF), 
         gnomad_minor = ifelse(AF > 0.5, REF, ALT),
         CHROM = as.numeric(CHROM)) %>%
  janitor::clean_names() %>%
  select(., chrom, pos, id, ref, alt, gnomad_minor, af, gnomad_maf) 

### snpEff 
system('java -Xmx8g -jar snpEff.jar GRCh37.75 examples/pdgwas_loci.vcf > output/pdgwas_loci.ann.vcf')

snpeff.path <- "intermediate/pdgwas_loci.ann.vcf"
snpeff.raw <- vcfR::read.vcfR(snpeff.path) %>%
  vcfR::vcfR2tidy(., info_only = F)

ann_cols <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO")
ann <- snpeff.raw$gt %>% 
  select(ChromKey, POS, SNP = gt_SNP) %>%
  left_join(., snpeff.raw$fix, by = c("ChromKey", "POS")) %>% 
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


## Combine 
nalls <- nalls_munged  %>%
  left_join(pd_gnomad_af, by = c("CHR" = 'chrom', 'BP' = 'pos', 'A1' = 'ref', 'A2' = 'alt')) %>%
  left_join(snpeff, by = c("SNP" = "snp")) 



nalls_raw %>% 
  select(1:10, failed_final_filtering_and_qc,locus_number) %>% 
  filter(failed_final_filtering_and_qc == 0) %>% 
  filter(qtl_nominated_gene_nearest_qtl %in% (nalls_smr_raw %>% filter(pass_bonferroni == "pass" & smr_p < 1e-8) %>% pull(gene)))

## Plots 

write_csv(nalls_munged, "intermediate/pd_loci.csv")

gene_labs = nalls_smr_raw %>% filter(pass_bonferroni == "pass" & smr_p < 1e-6) %>% pull(gene)
nalls.pdat <- nalls %>%
  mutate(
    dir = ifelse(OR > 1, "risk", "protective"), 
    effect = ifelse(OR > 1, OR, 1/OR),
    size = ifelse(effect < 1.2, 1, effect),
    architecture2 = case_when(
      MAF > 0.05 &  effect > 3 ~ "Common, High",
      MAF > 0.05 &  between(effect, 1.5, 3) ~ "Common, Intermediate",
      MAF > 0.05 &  between(effect, 1.1, 1.5) ~ "Common, Moderate",
      MAF > 0.05 &  effect < 1.1 ~ "Common, Low",
      between(MAF, 0.01, 0.05) &  effect > 3 ~ "Low, High",
      between(MAF, 0.01, 0.05) &  between(effect, 1.5, 3) ~ "Low, Intermediate",
      between(MAF, 0.01, 0.05) &  between(effect, 1.1, 1.5) ~ "Low, Moderate",
      between(MAF, 0.01, 0.05) &  effect < 1.1 ~ "Low, Low",
      between(MAF, 0.001, 0.01) &  effect > 3 ~ "Rare, High",
      between(MAF, 0.001, 0.01) &  between(effect, 1.5, 3) ~ "Rare, Intermediate",
      between(MAF, 0.001, 0.01) &  between(effect, 1.1, 1.5) ~ "Rare, Moderate",
      between(MAF, 0.001, 0.01) &  effect < 1.1 ~ "Rare, Low",
      MAF < 0.001 &  effect > 3 ~ "Very Rare, High",
      MAF < 0.001 &  between(effect, 1.5, 3) ~ "Very Rare, Intermediate",
      MAF < 0.001 &  between(effect, 1.1, 1.5) ~ "Very Rare, Moderate",
      MAF < 0.001 &  effect < 1.1 ~ "Very Rare, Low",
      TRUE ~ NA_character_
    ), 
   label = case_when( 
     SNP == "rs34637584" ~ "LRRK2 G2019S",
     SNP == "rs76763715" ~ "GBA1 N370S",
     SNP == "rs73038319" ~ QTL_NOMINATED_GENE_NEAREST_QTL,
     SNP == "rs117896735" ~ QTL_NOMINATED_GENE_NEAREST_QTL,
     QTL_NOMINATED_GENE_NEAREST_QTL %in% c("LRRK2") ~ QTL_NOMINATED_GENE_NEAREST_QTL,
     QTL_NOMINATED_GENE_NEAREST_QTL %in% gene_labs ~ QTL_NOMINATED_GENE_NEAREST_QTL,
     # !is.na(QTL_NOMINATED_GENE_NEAREST_QTL) ~ QTL_NOMINATED_GENE_NEAREST_QTL,
     TRUE ~ ""
   )
  )


# LRRK2, MAPT, and SNCA, GBA, 

theme.size = 8
geom.text.size = (theme.size - 2) * 0.36

### Bivariate Color scheme Legend
# out.legend <- cowplot::get_legend(gwas.p + theme(legend.title= element_blank()) + guides(color=guide_legend(ncol=3)))

bivariate_color_scale <- tribble(
  ~freq, ~effect, ~fill, ~architecture, 
  4, 4, "#453687", "Common, High", 
  4, 3, "#45589c", "Common, Intermediate",
  4, 2, "#477ab1", "Common, Moderate",
  4, 1, "#499cc6", "Common, Low",
  3, 4, "#713b8f", "Low, High",
  3, 3, "#7461a3", "Low, Intermediate",
  3, 2, "#7488ba", "Low, Moderate",
  3, 1, "#77add1", "Low, Low",
  2, 4, "#9e4096", "Rare, High",
  2, 3, "#9f6bad", "Rare, Intermediate",
  2, 2, "#a395c6", "Rare, Moderate",
  2, 1, "#a5bedc", "Rare, Low",
  1, 4, "#c9479d", "Very Rare, High",
  1, 3, "#cd74b4", "Very Rare, Intermediate",
  1, 2, "#cfa3ce", "Very Rare, Moderate",
  1, 1, "#d3cfe6", "Very Rare, Low"
)

use_col <- bivariate_color_scale$fill
names(use_col) <- bivariate_color_scale$architecture

bivariate_legend = bivariate_color_scale %>%
  ggplot(aes(x = freq, y = effect)) +
  geom_tile(aes(fill = fill)) +
  scale_fill_identity() + 
  scale_y_continuous(breaks = c(1,2,3,4), labels = c("Very\nLow", "Low", "Moderate", "High")) + 
  scale_x_continuous(breaks = c(1,2,3,4), labels = c("Very\nRare", "Rare", "Low", "Common")) + 
  theme_classic() + 
  labs(x = "Allele Frequency", y = "Effect Size") + 
  theme(
    text = element_text(size = 5),
    axis.line.y = element_line(
      arrow = grid::arrow(length = unit(0.1, "cm"), ends = "last")),
    axis.line.x = element_line(
      arrow = grid::arrow(length = unit(0.1, "cm"), ends = "last")),
    axis.text.x = element_text(),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  ) 

bivariate_legend

### Empty Y axis
pdyaxis.p <- ggplot() + 
  scale_x_continuous(limits = c(-0.1, 0.1), breaks=c(0), labels=c("Pathogenic Mutation")) +
  scale_y_continuous(trans='log', limits = c(0.1, 13.5), 
                     breaks=c(0.25, 0.5, 1, 2, 4, 8, 12), 
                     labels = c("0.25", "0.5", "1", "2", "4", "8", "12")) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  theme_classic() + 
  labs(x = "", y = "Odds Ratio - Minor Allele") + 
  theme(
    text = element_text(size = theme.size), 
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )

pdyaxis.p

# GWAS Architecture
pdgwas.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_point(data = nalls.pdat, aes(x = gnomad_maf, y = OR, color = architecture2, size = size)) +
  geom_text_repel(
    data = filter(nalls.pdat, label == "LRRK2 G2019S"),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, seed = 1,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(nalls.pdat, dir == "risk" & label != "LRRK2 G2019S"),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, nudge_y = 0.3, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(nalls.pdat, dir == "protective"),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, nudge_y = -0.3, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  scale_size(guide = 'none', range = c(0.5,6)) +
  theme_classic() + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log', limits = c(0.1, 13.5)) + 
  labs(x = "Global Population Minor Alelle Frequency", y = "Odds Ratio - Minor Allele") + 
  # scale_colour_manual(values = c("#3b4994", "#5ac8c8", "#8c62aa", "#dfb0d6", "#ace4e4", "red")) + 
  scale_colour_manual(values = use_col) + 
  guides(color=guide_legend(ncol=2)) + 
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.line.y=element_blank(), 
    text = element_text(size = theme.size), 
    legend.position = 'none', 
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )

pdgwas.p


### ADAD Genes
adpd <- tribble(
  ~Gene, ~x, ~y, 
  "SNCA", 1, 1,
  "PRKN", 1, 2,
  "PARK7", 1, 3,
  "LRRK2", 1, 4,
  "PINK1", 1, 5,
  "POLG", 2, 1.5,
  "ATP13A2", 2, 2.5,
  "FBXO7", 2, 3.5,
  "GBA", 2, 4.5,
  "PLA2G6", 3, 1,
  "VPS35", 3, 2,
  "DNAJC6", 3, 3,
  "SYNJ1", 3, 4,
  "VPS13C", 3, 5
)

pos <- ggbeeswarm::position_quasirandom()
adpd.p <- ggplot(adpd, aes(x = x, y = y, color = Gene, label = Gene)) + 
  geom_point(position = pos, size = 11, fill = "#be64ac", shape = 21, color = "#87497b") +
  geom_text(position = pos, color = "white", size = geom.text.size) +
  # geom_text_repel(data = filter(adpd, x == 1),
  #   aes(x = x, y = y, label = Gene, point.size = 5),
  #   max.overlaps = Inf, nudge_y = 0.15, seed = 333, xlim = c(1, 2),
  #   segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
  #   color = "black", size = geom.text.size, show.legend = F) +
  # geom_text_repel(data = filter(adpd, x == 2),
  #                 aes(x = x, y = y, label = Gene, point.size = 5),
  #                 max.overlaps = Inf, nudge_y = 0.15, seed = 333, xlim = c(2, 3),
  #                 segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
  #                 color = "black", size = geom.text.size, show.legend = F) +
  # geom_text_repel(data = filter(adpd, x == 3),
  #                 aes(x = x, y = y, label = Gene, point.size = 5),
  #                 max.overlaps = Inf, nudge_y = 0.15, seed = 333, xlim = c(3, 4),
  #                 segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
  #                 color = "black", size = geom.text.size, show.legend = F) +
  scale_y_continuous(breaks = 3, labels = "Causal") +
  theme_classic() +
  labs(x = "Minor Allele Frequency", y = " ") + 
  coord_cartesian(ylim=c(0, 5.5), xlim = c(0.3, 3.5)) +
  # geom_segment(x = -0.5, y = -0.5, xend = -0.5, yend = 1.5, color = 'black') +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.line.x=element_blank(), 
    axis.text.y = element_text(angle=90), 
    axis.ticks.y = element_blank(),
    # axis.line.y=element_blank(),
    text = element_text(size = theme.size)  
  )

adpd.p

## Plotgardner
png("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/plots/pd3_gwas.png", width = 9, height = 4.5, units = "in", res = 600)
tiff("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/plots/pd3_gwas.tiff", width = 9, height = 4.5, units = "in", res = 300)
pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = pdyaxis.p,
  x = 0, y = 0.75,
  width = 2, height = 3.75, just = c("left", "top")
)

plotGG(
  plot = adpd.p,
  x = 0.09, y = 0,
  width = 2, height = 1, just = c("left", "top")
)

plotGG(
  plot = pdgwas.p + theme(legend.position = "none"),
  x = 2, y = 0.75,
  width = 7, height = 3.75, just = c("left", "top")
)

plotGG(
  plot = bivariate_legend,
  x = 8.9, y = 0.1,
  width = 1.4, height = 1.4, just = c("right", "top")
)

pageGuideHide()
dev.off()











































