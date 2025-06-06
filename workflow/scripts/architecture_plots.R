# Genetic Architecture of Alzheimer's disease 
library(dplyr)
library(tidyr)
library(purrr)
library(janitor)
library(tibble)
library(glue)
library(stringr)
library(readr)
library(forcats)
library(glue)
library(ggrepel)
library(plotgardener)
`%nin%` = negate(`%in%`)

## Snakemake
### Input 
adgwas_loci.path = "results/adgwas_loci.csv"
adgwas_meta.path = "intermediate/gwas_metadata.csv"
adgwas_power.path = "intermediate/adgwas_power.csv"
bellenguez_supp.path = "resources/Bellenguez2022/41588_2022_1024_MOESM4_ESM.xlsx"

adgwas_loci.path = snakemake@input[['loci']]
adgwas_meta.path = snakemake@input[['meta']]
adgwas_power.path = snakemake@input[['power']]

### Output
outfile_abs.path = snakemake@output[['outpng_abs']]
outfile_obs.path = snakemake@output[['outpng_obs']]

## Import Dataset
message(
  "\nImporting: ",
  "\n\t", adgwas_loci.path,
  "\n\t", adgwas_meta.path,
  "\n\t", adgwas_power.path, 
  "\n"
)

adgwas_loci <- read_csv(adgwas_loci.path)
adgwas_meta <- read_csv(adgwas_meta.path)
adgwas_power <- read_csv(adgwas_power.path)

## Wrangle AD GWAS Loci datasets

### Gene labels to plot 
### Select variant with largest effective sample size in each locus 
### Remove APOE SNPS locus from AD datasets 
apoe_locus <- adgwas_loci %>% filter(SNP == "rs429358") %>% slice(1) %>% select(locus, cytoband, locus_ld)

adgwas_variants <- adgwas_loci %>% 
  # filter(study %in% c("Bellenguez", "Jonsson", "Reiman")) %>%
  arrange(locus_ld) %>% 
  mutate(locus_ld = ifelse(is.na(locus_ld), row(.), locus_ld)) %>%
  left_join(adgwas_meta, by = "study") %>%
  group_by(locus_ld) %>%
  arrange(-neff) %>%
  slice(1) %>%
  ungroup() %>% 
  filter(locus %nin% apoe_locus$locus) %>%
  bind_rows(filter(adgwas_loci, str_detect(SNP, "APOE")))

### Bellenguez Supplementary gene prioritization
st20 <- readxl::read_xlsx(bellenguez_supp.path, sheet = 21, skip = 2) %>% 
  janitor::clean_names() %>%
  select(gene, locus, gene_prioritization_tier, scoregene) %>%
  filter(gene_prioritization_tier %nin% c("-", "Excluded (Independent Loci)", "Tier 2", "Excluded (IGH Locus)")) %>%
  mutate(
    locus2 = str_extract(locus, ".*(?= Locus)")
  ) %>%
  select(-locus)

### Wrangling data for plotting 
### Lable rare variants with amino acid change, else selected gene lables
### Absolute effect sizes, direction, and architecture space
genes_to_plot <- c("ABCA7", "ABCA1", "TREM2", "APOE", "SIGLEC11", "GRN", "ABI3", "PLCG2", "BIN1", "CD2AP", "EED", "SPI1", 
                   "PTK2B", "CLU", "CASS4", "BLNK", "NCK2", "RIN3", "MS4A", "CD33", "PILRA", "CTSB", "CTSH", "APP", 
                   "CR1", "HLA", "EPHA1", "MS4A", "ADAM10", "MAPT", "ACE" 
                   )

dat.p <- adgwas_variants %>%
  left_join(st20, by = c("GENE" = "locus2")) %>%
  mutate(
    label = case_when(
      annotation_impact %in% c("HIGH", "LOW", "MODERATE") & gnomad_maf < 0.05 &  GENE %in% genes_to_plot ~ glue("{gene_name} {hgvs_p_new}"),
      str_detect(GENE, "APOE") ~ GENE, 
      GENE %in% genes_to_plot ~ GENE, 
      !is.na(gene) ~ GENE,
      # !is.na(gvc_gene) ~ gvc_gene,
    ), 
    label = ifelse(is.na(label), "", label),
    dir = ifelse(OR > 1, "risk", "protective"), 
    effect = ifelse(OR > 1, OR, 1/OR),
    size = ifelse(effect < 1.2, 1, effect),
    shp = ifelse(
      (label != "SORT1" & annotation_impact %in% c("HIGH", "LOW", "MODERATE") & 
        gnomad_maf < 0.05 & effect > 1.3 | str_detect(label, "APOE")), TRUE, FALSE
    ),
    category=cut(effect, breaks=c(1, 1.1, 1.2, 1.3, 1.5, 2, 2.5, 3, 3.5, Inf), 
                 labels=c("1-1.1","1.1-1.2","1.2-1.3","1.3-1.5","1.5-2", "2-2.5", "2.5-3", "3-3.5", ">4")),
    architecture = case_when(
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
      )
    ) %>% 
  select(SNP, CHR, BP, GENE, gnomad_maf, OR, label, shp, architecture, 
         annotation, annotation_impact, dir, size, effect, cytoband, locus, locus_ld, study) 

## Wrangle Power Curve datases 
power.dat <- adgwas_power %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or, 
         fill = "fill") %>% 
  # filter(study %in% c("Lambert", "Bellenguez")) %>%
  filter(study %in% c("Bellenguez", "Future1", "Future2")) %>%
  group_by(study) %>%
  fill(or, .direction = "up") %>%
  distinct(maf, or, .keep_all = T) %>%
  ungroup()
  
## Plots ======================================================================

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

ggsave("sandbox/plots/bivar_leg.png", plot = bivariate_legend + theme(legend.position = "none"), 
       width = 1.5, height = 1.5, units = "in")

### GWAS Architecture - observed scale
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
  geom_point(data = dat.p, aes(x = gnomad_maf, y = OR, color = architecture, size = size)) +
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

ggsave("sandbox/plots/AD_points.png", plot = adgwas.p + theme(legend.position = "none"), 
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
  ~Gene, ~x, ~y, 
  "APP", 1.5, 2,
  "PSEN1", 1, 1,
  "PSEN2", 2, 1,
)

pos <- ggbeeswarm::position_quasirandom(groupOnX = TRUE)
adad.p <- ggplot(adad, aes(x = x, y = y, color = Gene, label = Gene)) + 
  geom_point(position = pos, size = 11, fill = "#be64ac", shape = 21, color = "#87497b") +
  geom_text(position = pos, color = "white", size = geom.text.size) +
  scale_y_continuous(breaks = 1, labels = "Highly\nPenetrant") +
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
    axis.text.y = element_text(angle=90, hjust=0.5,vjust=0.5), 
    axis.ticks.y = element_blank(),
    # axis.line.y=element_blank(),
    text = element_text(size = theme.size), 
    panel.background = element_rect(fill='white'),
    plot.background = element_rect(fill='white', color=NA),
  )

adad.p

### Plotgardner
png(outfile_obs.path, width = 9, height = 4.5, units = "in", res = 600)
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

### Absolute Scale 
abs_gwas.p <- ggplot() + 
  # geom_line(data = power_0001_12.dat, aes(x = maf, y = or), lwd=0.25, col="#F66B0E") +
  geom_rect(data = tibble(ymax = 1.27, ymin = 0.98, xmax = 0.51, xmin = 0.049), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            color = "grey80", size = 0.25, fill = "grey95"
  ) + 
  geom_line(data = filter(power.dat, study == "Bellenguez"), aes(x = maf, y = or), lwd=0.25, color ="grey90") + 
  geom_line(data = filter(power.dat, study == "Future2"), aes(x = maf, y = or), lwd=0.25, color ="grey90") + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or) %>% rename(or = Bellenguez), 
              aes(ymin = or, ymax = Future2, x = maf), fill = "grey90") +
  # geom_hline(yintercept = 12) + 
  geom_text_repel(
    data =  dat.p %>% filter(gnomad_maf < 0.05),
    aes(x = gnomad_maf, y = effect, label = label, point.size = size),
    max.overlaps = Inf, seed = 100, xlim = c(NA, -1.30),
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F,) +
  geom_text_repel(
    data =  dat.p %>% filter(gnomad_maf > 0.05 & effect > 1.25),
    aes(x = gnomad_maf, y = effect, label = label, point.size = size),
    max.overlaps = Inf, seed = 100,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F,) +
  geom_point(data = filter(dat.p, dir == "risk" & shp == T) , shape = 24,
           aes(x = gnomad_maf, y = effect, 
               color = architecture, 
               fill = architecture, 
               size = size)) +
  geom_point(data = filter(dat.p, dir == "protective" & shp == T), shape = 25, 
             aes(x = gnomad_maf, y = effect, 
                 color = architecture, 
                 fill = architecture, 
                 size = size)) +
  geom_point(data = filter(dat.p, shp == F), shape = 19, 
             aes(x = gnomad_maf, y = effect, 
                 color = architecture, 
                 fill = architecture, 
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
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
    
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
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
  )
abs_yaxis.p


## The empty void 
void.p <- ggplot() + theme_void() + 
  theme(
    # panel.background = element_rect(fill='white'),
    plot.background = element_rect(fill='white', color=NA),
  )

png(outfile_abs.path, 
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



###############################################################################
# eBioMedicine

cutout.p <- ggplot() + 
  # geom_line(data = power_0001_12.dat, aes(x = maf, y = or), lwd=0.25, col="#F66B0E") +
  geom_smooth(data = filter(power.dat, study == "Bellenguez"), aes(x = maf, y = or), lwd=2.5, color ="grey90", se = F) + 
  geom_smooth(data = filter(power.dat, study == "Future2"), aes(x = maf, y = or), lwd=2.5, color ="grey90", se = F) + 
  geom_ribbon(data = select(power.dat, maf, or, study) %>% pivot_wider(names_from = study, values_from = or) %>% rename(or = Bellenguez), 
              aes(ymin = or, ymax = Future2, x = maf), fill = "grey90") +
  geom_text_repel(
    data = dat.p %>% filter(gnomad_maf > 0.05 & effect < 1.25),
    aes(x = gnomad_maf, y = effect, label = label, point.size = size),
    max.overlaps = Inf, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    angle = 90, ylim = c(NA, 0.99), hjust = 1,
    color = "black", size = geom.text.size, show.legend = F, 
    direction = "x",
    segment.curvature = -1e-20, segment.ncp = 3) +
  geom_point(data = dat.p,
             aes(x = gnomad_maf, y = effect, 
                 color = architecture, 
                 fill = architecture), 
             size = 0.5) +
  scale_size(guide = 'none', range = c(0.5,6)) +
  theme_light() + 
  scale_x_continuous(#trans='log',
    breaks = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("0.05", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(breaks = c(1, 1.1, 1.2)) +
  coord_cartesian(xlim=c(0.049, 0.5), ylim = c(0.85, 1.2)) + 
  # scale_colour_manual(values = c("#3b4994", "#5ac8c8", "#8c62aa", "#dfb0d6", "#ace4e4", "red")) + 
  scale_colour_manual(values = use_col) + 
  scale_fill_manual(values = use_col) + 
  guides(color=guide_legend(ncol=2)) + 
  theme(
    axis.title.x = element_blank(),
    # axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.line.y=element_blank(), 
    text = element_text(size = theme.size), 
    legend.position = 'none', 
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
  )

cutout.p

## Zoom-in polygon 
positions = tribble(
  ~x, ~y, 
  0.1, 1.34,
  7, 2,
  # 6.75, 1.85,
  # 8.65, 1.85, 
  9.05, 2, 
  9.36, 1.34
)

ribbon.p <- ggplot() +
  geom_polygon(data = positions, aes(x = x, y = y), fill = 'grey95') + 
  theme_void() +
  coord_cartesian(xlim=c(0, 9), ylim = c(0, 4.5)) + 
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
  )


tiff(
    filename = "results/plots/eBioMedicine_AD_GeneticArchitecture.tiff",
    width    = 9,           # inches
    height   = 4.5,         # inches
    units    = "in",
    res      = 300,          # will embed 300×300 dpi
    compression = "lzw", type = "cairo"
  )

pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = ribbon.p,
  x = 0, y = 0,
  width = 9, height = 4.5, just = c("left", "top")
)

plotGG(
  plot = abs_gwas.p + theme(legend.position = "none"),
  x = 0.1, y = 0.3,
  width = 9, height = 2.7, just = c("left", "top")
)

plotGG(
  plot = cutout.p + theme(legend.position = "none"),
  x = 0.24, y = 3,
  width =  8.76, height = 1.5, just = c("left", "top")
)

plotGG(
  plot = void.p,
  x = 0, y = 0.75,
  width = 1.5, height = 2, just = c("left", "top")
)

plotGG(
  plot = abs_yaxis.p,
  x = 0, y = 0.75,
  width = 1.5, height = 2.25, just = c("left", "top")
)

# plotGG(
#   plot = void.p,
#   x = 0, y = 0.75,
#   width = 0.6, height = 0.25, just = c("left", "top")
# )

plotGG(
  plot = adad.p,
  x = 0.04, y = 0,
  width = 1.5, height = 0.75, just = c("left", "top")
)

# plotGG(
#   plot = abs_gwas.p + theme(legend.position = "none"),
#   x = 1.5, y = 0.75,
#   width = 7.5, height = 3.75, just = c("left", "top")
# )

plotGG(
  plot = bivariate_legend,
  x = 8.9, y = 0,
  width = 1.4, height = 1.4, just = c("right", "top")
)

pageGuideHide()
dev.off()

library(magick)

file_in  <- "results/plots/eBioMedicine_AD_GeneticArchitecture.tiff"
img <- image_read(file_in)
image_info(img)






























