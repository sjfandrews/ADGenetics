library(MungeSumstats)
library(tidyverse)
library(tabulizer)
library(fuzzyjoin)
library(plotgardener)

data(sumstatsColHeaders)
setwd("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/")

chia_ss.raw <- read_tsv("~/Downloads/GCST90001390_buildGRCh38.tsv.gz")

# ~/Dropbox/Research/PostDoc-MSSM/ADGenetics/resources/Chia2021.pdf

chia_raw <- chia_ss.raw %>%
  filter(variant_id %in% c("rs2230288", "rs6733839", "rs6599388", "rs7680557", "rs769449"))
  
  
chia_clean <- chia_raw %>%
  rename(snp = variant_id, chr = chromosome, bp = base_pair_location, or = odds_ratio) %>%
  MungeSumstats::format_sumstats(path=., 
                                 ref_genome="GRCh38", 
                                 convert_ref_genome = "GRCh37",
                                 allele_flip_check = FALSE, 
                                 allele_flip_drop = FALSE, 
                                 bi_allelic_filter = FALSE,
                                 return_data = TRUE, 
                                 log_folder = "data/MungeSumstats", 
                                 force_new = TRUE) %>%
  mutate(
    MAF = ifelse(FRQ > 0.5, 1 - FRQ, FRQ), 
    BETA = ifelse(FRQ > 0.5, BETA * -1, BETA), 
    OR = ifelse(FRQ > 0.5, 1/OR, OR)
  ) 
  
write_csv(chia_clean, "intermediate/lbd_loci.csv")

chia_clean %>%
  select(-SNP) %>%
  MungeSumstats::write_sumstats(., "intermediate/lbdgwas_loci.vcf", write_vcf =T)

## Gnomad MAF
system('bcftools view -R input/lbdgwas_loci.vcf -Ov -o output/lbdgwas_loci_maf.vcf reference/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz')

gnomadaf.raw <- vcfR::read.vcfR('~/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/lbdgwas_loci_maf.vcf') 

gnomad_af <- gnomadaf.raw %>%
  vcfR::vcfR2tidy(., info_only = T) %>%
  magrittr::extract2("fix") %>%
  select(CHROM, POS, ID, REF, ALT, AF, AF_raw, AF_afr, AF_eas, AF_amr, AF_nfe, AF_asj, AF_fin) %>%
  mutate(AF = as.numeric(AF), 
         gnomad_maf = ifelse(AF > 0.5, 1 - AF, AF), 
         gnomad_minor = ifelse(AF > 0.5, REF, ALT),
         CHROM = as.character(CHROM)) %>%
  janitor::clean_names() %>%
  select(., chrom, pos, id, ref, alt, gnomad_minor, af, gnomad_maf) %>%
  slice(-5) 


## Combind data 
lbd.pdat <- chia_clean %>%
  left_join(select(gnomad_af, chrom, pos, af, gnomad_maf), by = c("CHR" = 'chrom', 'BP' = 'pos')) %>%
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
      SNP == "rs2230288" ~ "GBA E326K",
      SNP == "rs6733839" ~ "BIN1",
      SNP == "rs6599388" ~ "TMEM175",
      SNP == "rs7680557" ~ "SNCA-AS1",
      SNP == "rs769449" ~ "APOE ε4",
      TRUE ~ ""
    )
  )
  


## Plots

# GWAS Plot
lbd.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_point(data = lbd.pdat, aes(x = gnomad_maf, y = OR, color = architecture2, size = size)) +
  geom_text_repel(
    data = lbd.pdat %>% filter(SNP %nin% c("rs2230288", "rs769449")),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = lbd.pdat %>% filter(SNP %in% c("rs2230288", "rs769449")),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, nudge_y = 0.25, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  scale_size(guide = 'none', range = c(0.5,6)) +
  theme_classic() + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log', limits = c(0.1, 13.5)) + 
  labs(x = "Global Population Minor Allele Frequency", y = "Odds Ratio - Minor Allele") + 
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

lbd.p

## Dominant mutations
adlbd <- tribble(
  ~Gene, ~x, ~y, 
  "SNCA", 1, 1,
  "GBA", 2, 2,
)

pos <- ggbeeswarm::position_quasirandom()
adlbd.p <- ggplot(adlbd, aes(x = x, y = y, color = Gene, label = Gene)) + 
  geom_point(position = pos, size = 11, fill = "#be64ac", shape = 21, color = "#87497b") +
  geom_text(position = pos, color = "white", size = geom.text.size) +
  scale_y_continuous(breaks = 2, labels = "Causal") +
  theme_classic() +
  labs(x = "Global Population Minor Allele Frequency", y = " ") + 
  coord_cartesian(ylim=c(0, 3), xlim = c(0.5, 2.5)) +
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

adlbd.p

## y axis
lbdyaxis.p <- ggplot() + 
  scale_x_continuous(limits = c(-0.1, 0.1), breaks=c(0), labels=c("Pathogenic Mutation")) +
  scale_y_continuous(trans='log', limits = c(0.1, 13), 
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

png("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/plots/lbd_gwas.png", width = 9, height = 4.5, units = "in", res = 600)
tiff("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/plots/lbd_gwas.tiff", width = 9, height = 4.5, units = "in", res = 300)
pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = lbdyaxis.p,
  x = 0, y = 0.75,
  width = 2, height = 3.75, just = c("left", "top")
)

plotGG(
  plot = adlbd.p,
  x = 0.09, y = 0,
  width = 2, height = 1, just = c("left", "top")
)

plotGG(
  plot = lbd.p + theme(legend.position = "none"),
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
















