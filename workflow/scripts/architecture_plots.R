library(tidyverse)
library(glue)
library(ggrepel)
library(plotgardener)
`%nin%` = negate(`%in%`)

setwd("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADgenetics/")

ad_loci <- read_csv("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/ad_loci.csv")
adgwas <- read_csv("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/ad_gwas_meta_data.csv")

## AD GWAS Loci Plots

genes_to_plot <- c("ABCA7", "ABCA1", "TREM2", "APOE", "SIGLEC11", "GRN", "ABI3", "PLCG2", "BIN1", "CD2AP", "EED", "SPI1", 
"PTK2B", "CLU", "CASS4", "BLNK", "NCK2", "RIN3", "MS4A", "CD33", "PILRA", "CTSB", "CTSH", "APP")

ad_loci %>% 
  left_join(adgwas, by = "study") %>%
  group_by(locus) %>%
  arrange(-neff) %>%
  slice(1) %>%
  ungroup() %>%
  filter(., annotation_impact %in% c("HIGH", "LOW", "MODERATE")) %>% 
  filter(gnomad_maf < 0.05) %>% 
  select(locus, cytoband, SNP, CHR, BP, Major, Minor, allele, gnomad_maf, OR, P, gencode_gene, gvc_gene, gene_name, annotation, annotation_impact, hgvs_p) 
  
apoe_locus <- ad_loci %>% filter(SNP == "rs429358") %>% slice(1) %>% select(locus, locus_ld)
  
dat_loci <- ad_loci %>% 
  arrange(locus_ld) %>% 
  mutate(locus_ld = ifelse(is.na(locus_ld), row(.), locus_ld)) %>%
  left_join(adgwas, by = "study") %>%
  group_by(locus_ld) %>%
  arrange(-neff) %>%
  slice(1) %>%
  ungroup() %>% 
  filter(locus %nin% "chr19:44851516-46741841") %>%
  mutate( # TOFIX 
    OR = ifelse(SNP == "rs139643391", exp(-0.0619), OR),
    gnomad_maf = ifelse(SNP == "rs149080927", 0.3829, gnomad_maf),
    gnomad_maf = ifelse(str_detect(SNP, "APOE"), FRQ, gnomad_maf),
    locus = ifelse(str_detect(SNP, "APOE"), apoe_locus$locus, locus),
    locus_ld = ifelse(str_detect(SNP, "APOE"), apoe_locus$locus_ld, locus_ld)
    ) %>%
  mutate(
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
  mutate(aa = str_extract(hgvs_p,"\\d{1,4}"), 
         hgvs_p_new = case_when(
           Minor == allele ~ glue("{p1}{aa}{p2}"),
           Minor != allele ~ glue("{p2}{aa}{p1}")
         )
  ) %>%
  mutate(
    label = case_when(
      annotation_impact %in% c("HIGH", "LOW", "MODERATE") & gnomad_maf < 0.05 &  GENE %in% genes_to_plot ~ glue("{gene_name} {hgvs_p_new}"),
      str_detect(GENE, "APOE") ~ GENE, 
      GENE %in% genes_to_plot ~ GENE
      # !is.na(gvc_gene) ~ gvc_gene,
    ),
    OR = ifelse(SNP == "rs429358", 3.6, OR),
    dir = ifelse(OR > 1, "risk", "protective"), 
    effect = ifelse(OR > 1, OR, 1/OR),
    category=cut(effect, breaks=c(1, 1.1, 1.2, 1.3, 1.5, 2, 2.5, 3, 3.5, Inf), 
                 labels=c("1-1.1","1.1-1.2","1.2-1.3","1.3-1.5","1.5-2", "2-2.5", "2.5-3", "3-3.5", ">4")),
    architecture = case_when(
      (OR >= 1.6 | OR <= 0.6) & gnomad_maf < 0.001 ~ "Rare, Moderate",
      (OR >= 1.6 | OR <= 0.6) & gnomad_maf < 0.05 ~ "Low, Moderate",
      (OR < 1.6 & OR > 0.6) & gnomad_maf < 0.05 ~ "Low, Low",
      (OR > 1.6 | OR < 0.6) & gnomad_maf > 0.05 ~ "Common, High",
      (OR < 1.6 & OR > 0.6) & gnomad_maf > 0.05 ~ "Common, Low",
      TRUE ~ NA_character_
    ), 
    architecture2 = case_when(
      gnomad_maf > 0.05 &  effect > 3 ~ "Common, High",
      gnomad_maf > 0.05 &  between(effect, 1.5, 3) ~ "Common, Intermediate",
      gnomad_maf > 0.05 &  between(effect, 1.1, 1.5) ~ "Common, Moderate",
      gnomad_maf > 0.05 &  effect < 1.1 ~ "Common, Low",
      between(gnomad_maf, 0.01, 0.05) &  effect > 3 ~ "Low, High",
      between(gnomad_maf, 0.01, 0.05) &  between(effect, 1.5, 3) ~ "Low, Intermediate",
      between(gnomad_maf, 0.01, 0.05) &  between(effect, 1.1, 1.5) ~ "Low, Moderate",
      between(gnomad_maf, 0.01, 0.05) &  effect < 1.1 ~ "Low, Low",
      between(gnomad_maf, 0.001, 0.01) &  effect > 3 ~ "Rare, High",
      between(gnomad_maf, 0.001, 0.01) &  between(effect, 1.5, 3) ~ "Rare, Intermediate",
      between(gnomad_maf, 0.001, 0.01) &  between(effect, 1.1, 1.5) ~ "Rare, Moderate",
      between(gnomad_maf, 0.001, 0.01) &  effect < 1.1 ~ "Rare, Low",
      gnomad_maf < 0.001 &  effect > 3 ~ "Very Rare, High",
      gnomad_maf < 0.001 &  between(effect, 1.5, 3) ~ "Very Rare, Intermediate",
      gnomad_maf < 0.001 &  between(effect, 1.1, 1.5) ~ "Very Rare, Moderate",
      gnomad_maf < 0.001 &  effect < 1.1 ~ "Very Rare, Low",
      TRUE ~ NA_character_
      ),
    snpeff = case_when(
      (OR < 1.3 & OR > 0.75) & gnomad_maf > 0.05 ~ "Low effect size common variant",
      annotation == "downstream_gene_variant" ~ 'intergenic_region',
      annotation == "upstream_gene_variant" ~ 'intergenic_region',
      annotation == "regulatory_region_variant" ~ 'intergenic_region',
      annotation == "intron_variant,non_coding_transcript_variant" ~ 'intron_variant',
      annotation == "3_prime_UTR_variant" ~ 'UTR_variant',
      annotation == "5_prime_UTR_variant" ~ 'UTR_variant',
      annotation == "stop_gained" ~ 'missense_variant',
      annotation == "conservative_inframe_deletion" ~ 'missense_variant',
      TRUE ~ annotation
    ),
  ) 
 
filter(dat.p, architecture != "Other") %>% select(SNP, GENE, OR, FRQ, architecture)
filter(dat.p, str_detect(GENE, "APOE")) %>% select(SNP, OR, FRQ, architecture)

filter(dat.p, (OR >= 1.3 | OR <= 0.75) & FRQ <= 0.05)
filter(dat.p, (OR < 1.3 & OR > 0.75) & FRQ < 0.05)
filter(dat.p, between(OR, 1.3, 0.75) & FRQ < 0.05)
filter(dat.p, (OR > 1.3 | OR < 0.75) & FRQ > 0.05)
filter(dat.p, (OR < 1.3 & OR > 0.75) & FRQ > 0.05)

filter(dat.p, GENE == "SORL1")
filter(dat.p, str_detect(GENE, "TMEM106B"))
filter(dat.p, architecture == "Low frequency, Small effect") %>%
  select(SNP, GENE, FRQ, OR, study)

dat_loci %>% 
  select(locus, cytoband, locus_ld, SNP, CHR, BP, A1, A2, GENE, gnomad_maf, OR, label, architecture, architecture2, annotation, annotation_impact, dir, effect, study, gencode_gene, gencode_dist, gvc_gene, gvc_dist) %>%
  filter(!is.na(gvc_gene)) %>%
  group_split(gvc_gene)

################# Draw curves indicating the effect size needed across different MAF threshols
power.raw <- read_csv("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/ad_power_matti.csv") 

power.dat <- power.raw %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or, 
         fill = "fill") %>% 
  # filter(study %in% c("Lambert", "Bellenguez")) %>%
  filter(study %in% c("Bellenguez", "Future1", "Future2")) %>%
  group_by(study) %>%
  fill(or, .direction = "up") %>%
  distinct(maf, or, .keep_all = T) %>%
  ungroup()
  
 
test <- select(power.dat, maf, or, study) %>% 
  pivot_wider(names_from = study, values_from = or) 

ggplot(test, aes(y = Bellenguez, x = maf)) + 
  geom_line(data = test, aes(y = Bellenguez), lwd=0.25, color ="grey70") + 
  geom_line(data = test, aes(y = Lambert), lwd=0.25, color ="grey70") + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or), 
              aes(ymin = Bellenguez, ymax = Lambert), fill = "grey70") +
  scale_x_continuous(trans='log10', breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.00001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  theme_light()
  
  
ggplot(power.dat, aes(y = or, x = maf)) + 
  geom_line(data = filter(power.dat, study == "Bellenguez"), lwd=0.25, color ="grey70") + 
  geom_line(data = filter(power.dat, study == "Lambert"), lwd=0.25, color ="grey70") + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or) %>% rename(or = Bellenguez), 
              aes(ymin = or, ymax = Lambert), fill = "grey70") +
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  theme_light()

ggplot(power.dat, aes(y = or, x = maf)) + 
  geom_line(data = filter(power.dat, study == "Bellenguez"), lwd=0.25, color ="grey70") + 
  # geom_line(data = filter(power.dat, study == "Future1"), lwd=0.25, color ="grey70") + 
  geom_line(data = filter(power.dat, study == "Future2"), lwd=0.25, color ="grey70") + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or) %>% rename(or = Bellenguez),
              aes(ymin = or, ymax = Future2, alpha = maf), fill = "grey70") +
  scale_x_continuous(trans='log10', breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.00001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  # scale_alpha()
  theme_light()


## Plots ======================================================================

dat.p <- dat_loci %>% 
  # filter(study %in% c("Bellenguez", "Jonsson", "Reiman")) %>%
  select(SNP, CHR, BP, GENE, gnomad_maf, OR, label, architecture, architecture2, 
         annotation, annotation_impact, dir, effect, cytoband, locus, locus_ld, study) %>%
  mutate(size = ifelse(effect < 1.2, 1, effect), 
         label = ifelse(is.na(label), "", label)) 

theme.size = 8
geom.text.size = (theme.size - 2) * 0.36

### GWAS Plot 

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
  scale_y_continuous(breaks = c(1,2,3,4), labels = c("Low", "Modest", "Intermt.", "High")) + 
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
ggsave("~/Downloads/bivar_leg.png", plot = bivariate_legend + theme(legend.position = "none"), 
       width = 1.5, height = 1.5, units = "in")

ggplot(power.dat, aes(y = or, x = maf)) + 
  geom_line(data = filter(power.dat, study == "Bellenguez"), lwd=0.25, color ="grey70") + 
  geom_line(data = filter(power.dat, study == "Lambert"), lwd=0.25, color ="grey70") + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or) %>% rename(or = Bellenguez), 
              aes(ymin = or, ymax = Lambert), fill = "grey70") +
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  theme_light()

# GWAS Architecture
adgwas.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  # geom_line(data = power_0001_12.dat, aes(x = maf, y = or), lwd=0.25, col="#F66B0E") +
  # geom_line(data = power_0001_12.dat, aes(x = maf, y = inv_or), lwd=0.25, col="#F66B0E") +
  geom_text_repel(
    data = filter(dat.p, dir == "risk"),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, nudge_y = 0.15, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(dat.p, dir == "protective"),
    aes(x = gnomad_maf, y = OR, label = label, point.size = size),
    max.overlaps = Inf, nudge_y = -0.25, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_point(data = dat.p, aes(x = gnomad_maf, y = OR, color = architecture2, size = size)) +
  scale_size(guide = 'none', range = c(0.5,6)) +
  theme_classic() + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log', limits = c(0.1, 13)) + 
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
    legend.position = 'none'
  )

adgwas.p

ggsave("~/Downloads/AD_points.png", plot = gwas.p + theme(legend.position = "none"), 
       width = 7.5, height = 3.75, units = "in")

### Empty Y axis
adyaxis.p <- ggplot() + 
  scale_x_continuous(limits = c(-0.1, 0.1), breaks=c(0), labels=c("Pathogenic Mutation")) +
  scale_y_continuous(trans='log', limits = c(0.1, 13), 
                     breaks=c(0.25, 0.5, 1, 2, 4, 8, 12), 
                     labels = c("0.25", "0.5", "1", "2", "4", "8", "12")) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  theme_classic() + 
  labs(x = "", y = "Odds Ratio - Minor Allele") + 
  theme(
    text = element_text(size = theme.size)
  )

adyaxis.p

### ADAD Genes
adad <- tribble(
  ~Gene, ~OR, ~FRQ, 
  "APP", 1.2, 0.5,
  "PSEN1", 0.2, 0,
  "PSEN2", 0, 1,
)


adad <- tribble(
  ~Gene, ~x, ~y, 
  "APP", 1.5, 2,
  "PSEN1", 1, 1,
  "PSEN2", 2, 1,
)

pos <- ggbeeswarm::position_quasirandom()
adad.p <- ggplot(adad, aes(x = x, y = y, color = Gene, label = Gene)) + 
  geom_point(position = pos, size = 11, fill = "#be64ac", shape = 21, color = "#87497b") +
  geom_text(position = pos, color = "white", size = geom.text.size) +
  scale_y_continuous(breaks = 2, labels = "Causal") +
  theme_classic() +
  labs(x = "Global Population Minor Allele Frequency", y = " ") + 
  coord_cartesian(ylim=c(-0.5, 3), xlim = c(0.5, 2.5)) +
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
    text = element_text(size = theme.size), 
    panel.background = element_rect(fill='white'),
    plot.background = element_rect(fill='white', color=NA),
    
  )

adad.p

## Plotgardner
png("~/Dropbox/Research/PostDoc-MSSM/ADgenetics/plots/AD_GWAS2_BellenguezOnly.png", width = 9, height = 4.5, units = "in", res = 600)
tiff("~/Dropbox/Research/PostDoc-MSSM/ADgenetics/plots/AD_GWAS2_BellenguezOnly.tiff", width = 9, height = 4.5, units = "in", res = 300)

pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = adad.p,
  x = 0.09, y = 0,
  width = 2, height = 1, just = c("left", "top")
)

plotGG(
  plot = adyaxis.p,
  x = 0, y = 0.75,
  width = 2, height = 3.75, just = c("left", "top")
)

plotGG(
  plot = adgwas.p + theme(legend.position = "none"),
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

## Cowplot
cowplot::plot_grid(adad.p, gg, yaxis.p, gwas.p + theme(legend.position = "none"), 
                   # axis = "lb", align	= "hv",
                   nrow = 2, ncol = 2, 
                   rel_heights = c(0.2, 1),
                   rel_widths = c(0.2, 1))
ggsave("~/Downloads/AD_GWAS.png", width = 9, height = 4.5, units = "in")


### Absolute Scale 

abs_gwas.p <- ggplot() + 
  # geom_line(data = power_0001_12.dat, aes(x = maf, y = or), lwd=0.25, col="#F66B0E") +
  geom_line(data = filter(power.dat, study == "Bellenguez"), aes(x = maf, y = or), lwd=0.25, color ="grey90") + 
  geom_line(data = filter(power.dat, study == "Future2"), aes(x = maf, y = or), lwd=0.25, color ="grey90") + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or) %>% rename(or = Bellenguez), 
              aes(ymin = or, ymax = Future2, x = maf), fill = "grey90") +
  # geom_hline(yintercept = 12) + 
  geom_text_repel(
    data = dat.p,
    aes(x = gnomad_maf, y = effect, label = label, point.size = size),
    max.overlaps = Inf, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_point(data = filter(dat.p, dir == "risk" & label != "") , shape = 24,
           aes(x = gnomad_maf, y = effect, 
               color = architecture2, 
               fill = architecture2, 
               size = size)) +
  geom_point(data = filter(dat.p, dir == "protective" & label != ""), shape = 25, 
             aes(x = gnomad_maf, y = effect, 
                 color = architecture2, 
                 fill = architecture2, 
                 size = size)) +
  geom_point(data = filter(dat.p, label == ""), shape = 19, 
             aes(x = gnomad_maf, y = effect, 
                 color = architecture2, 
                 fill = architecture2, 
                 size = size)) +
  scale_size(guide = 'none', range = c(0.5,6)) +
  theme_classic() + 
  scale_x_continuous(trans='log10',
                     breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("5e-1", "0.0001", "0.001", "0.01", "0.1", "0.25", "0.5")) + 
  scale_y_continuous(trans='log') +
  labs(x = "Global Population Minor Alelle Frequency", y = "Odds Ratio - Minor Allele") + 
  coord_cartesian(xlim=c(0.00002, 0.5), ylim = c(0.9, 25)) + 
  # scale_colour_manual(values = c("#3b4994", "#5ac8c8", "#8c62aa", "#dfb0d6", "#ace4e4", "red")) + 
  scale_colour_manual(values = use_col) + 
  scale_fill_manual(values = use_col) + 
  guides(color=guide_legend(ncol=2)) + 
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.line.y=element_blank(), 
    text = element_text(size = theme.size), 
    legend.position = 'none', 
  )

abs_gwas.p

### Empty Y axis
abs_yaxis.p <- ggplot() + 
  # geom_hline(yintercept = 12) + 
  scale_x_continuous(limits = c(-0.1, 0.1), breaks=c(0), labels=c("Pathogenic Mutation")) +
  scale_y_continuous(trans='log', limits = c(0.9, 13), 
                     breaks=c(0, 1, 2, 4, 8, 12), 
                     labels = c("0", "1", "2", "4", "8", "12")) +
  # geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  theme_classic() + 
  labs(x = "", y = "Odds Ratio - Minor Allele\n(Absolute Scale)") + 
  theme(
    text = element_text(size = theme.size), 
    panel.background = element_rect(fill='white'),
    plot.background = element_rect(fill='white', color=NA),
    
  )
abs_yaxis.p


## The empty void 
void.p <- ggplot() + theme_void() + 
  theme(
    # panel.background = element_rect(fill='white'),
    plot.background = element_rect(fill='white', color=NA),
  )

png("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/plots/AD_GWASabs2_area.png", 
    width = 9, height = 4.5, units = "in", res = 300)
pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = abs_gwas.p + theme(legend.position = "none"),
  x = 0.1, y = 0,
  width = 9, height = 4.5, just = c("left", "top")
)

plotGG(
  plot = void.p,
  x = 0, y = 4,
  width = 1.5, height = 0.5, just = c("left", "top")
)

plotGG(
  plot = abs_yaxis.p,
  x = 0, y = 0.75,
  width = 1.5, height = 3.75, just = c("left", "top")
)

plotGG(
  plot = void.p,
  x = 0, y = 0.75,
  width = 0.6, height = 0.25, just = c("left", "top")
)

plotGG(
  plot = adad.p,
  x = 0.14, y = 0,
  width = 1.5, height = 1, just = c("left", "top")
)

# plotGG(
#   plot = abs_gwas.p + theme(legend.position = "none"),
#   x = 1.5, y = 0.75,
#   width = 7.5, height = 3.75, just = c("left", "top")
# )

plotGG(
  plot = bivariate_legend,
  x = 8.9, y = 0.1,
  width = 1.4, height = 1.4, just = c("right", "top")
)


pageGuideHide()
dev.off()































