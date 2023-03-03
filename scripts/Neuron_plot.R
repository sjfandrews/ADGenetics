## AD Genetic Architecture for Neuron review

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
# setwd("~/Dropbox/Research/PostDoc-MSSM/ADGenetics")
adgwas_loci.path = "results/adgwas_loci.csv"
adgwas_meta.path = "intermediate/gwas_metadata.csv"
adgwas_power.path = "intermediate/adgwas_power.csv"

adgwas_loci.path = snakemake@input[['loci']]
adgwas_meta.path = snakemake@input[['meta']]
adgwas_power.path = snakemake@input[['power']]

### Output
outfile_neuron.path = '~/Dropbox/Research/PostDoc-MSSM/ADGenetics/results/plots/Neuron_AD_GeneticArchitecture.png'
outfile_neuron.path = snakemake@output[['outpng_abs']]

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

## Wrangle data 
### AD GWAS meta data 
adgwas <- adgwas_meta %>% 
  filter(study %in% c("Lambert", "Kunkle", "Jansen", "Wightman", "Bellenguez")) %>%
  mutate(
    study = fct_relevel(study, c("Lambert", "Jansen", "Kunkle", "Wightman", "Bellenguez")), 
  ) %>%
  select(study, year, n_loci, neff, n) %>%
  pivot_longer(c("neff", "n"), names_to = "model", values_to = "size") 

## From Bellenguez
microglia_gene_list <- c("SPI1", "MAF", "MEF2C", "NR1H2", "NR1H3", "TREM2", "BLNK", "MS4A4A", "MS4A6A", 
                         "CD33", "SIGLEC11", "INPP5D/SHIP1", "NCK2", "ABI3", "PICALM", "EED", "CD2AP", 
                         "BIN1", "SORL1", "RAB10", "SPPL2A", "PGRN", "ZYX", "RIN3", "AP4E1", "AP4M1", 
                         "TP3INP1", "APOE", "CLU", "LACTB", "OPA3")

st20 <- readxl::read_xlsx("resources/Bellenguez2022/41588_2022_1024_MOESM4_ESM.xlsx", sheet = 21, skip = 2) %>% 
  janitor::clean_names() %>%
  select(gene, locus, ends_with("percent")) %>%
  rename(locus_bellenguez = locus) %>%
  mutate_at(vars(ends_with("percent")), as.numeric)

# st20 <- readxl::read_xlsx("resources/Bellenguez2022/41588_2022_1024_MOESM4_ESM.xlsx", sheet = 21, skip = 2) %>% 
#   janitor::clean_names() %>%
#   select(gene, gene_prioritization_tier, scoregene, scoretop) 

## Gene label appears in text
### R0
# intext_list <- c("IDUA", "JAZF1", "TSPAN14", "BLNK", "PLEKHA1", "CTSH", "SIGLEC11", 
#             "CR1", "BIN1", "INPP5D", "HLA-DQA1", "TREM2", "EPDR1", "PTK2B", 
#             "USP6NL", "EED", "PLCG2", "ABI3", "TSPOAP1", "SPPL2A", "SPI1", "ABCA1", 
#             "MS4A4A", "SORL1", "MAF"
# )

### R1
intext_list <- c('TREM2', 'PLCG2', 'BLNK', 'MS4A4A', 'MS4A6A', 'PILRA', 'SIGLEC11', 'APOE', 'CLU', 'NCK2', 'EED',
                 'ABI3', 'IL34', 'PICALM', 'BIN1', 'SORL1', 'RAB10', 'SPPL2A', 'GRN', 'ZYX', 'RIN3', 'AP4E1',
                 'AP4M1', 'TP53INP1', 'CELF1', 'SPI1', 'MAF', 'NR1H3', 'NR1H2', 'ABCA1', 'ABCA7', 'LACTB', 'OPA3', 
                 'EIF2B2', 'INPPFD', 'RHOH', 'SNX1', 'PLIN2'
)



## Microglia signature gene lists 
gosselin.raw <- readxl::read_xlsx("resources/gosselin/aal3222_gosselin_tables2.xlsx", sheet = 2) %>% 
  janitor::clean_names()

gosselin <- gosselin.raw %>% 
  mutate(microglia_gene = TRUE)

## Wrangle AD Loci
dat.p <- adgwas_loci %>% 
  filter(study %in% c("Bellenguez", "Reiman")) %>%
  filter(SNP != "rs429358") %>%
  left_join(gosselin, by = c("GENE" = "gene_name")) %>%
  left_join(st20, by = c("GENE" = "gene")) %>%
  mutate(
    dir = ifelse(OR > 1, "risk", "protective"), 
    effect = ifelse(OR > 1, OR, 1/OR),
    size = ifelse(effect < 1.2, 1, effect), 
       microglia_gene = ifelse(is.na(microglia_gene), FALSE, microglia_gene),
       intext = ifelse(GENE %in% intext_list | str_detect(GENE, "APOE"), TRUE, FALSE), 
       # microglia_gene = case_when(
       #   GENE %in% microglia_gene_list ~ TRUE,
       #   TRUE ~ FALSE),
       label = case_when(
         annotation_impact %in% c("HIGH", "LOW", "MODERATE") & gnomad_maf < 0.05  & microglia_gene == TRUE ~ glue("{gene_name} {hgvs_p_new}"),
         str_detect(GENE, "APOE") ~ GENE,
         GENE == 'SLC24A4' ~ "RIN3", 
         GENE == 'MS4A4A' ~ 'MS4A4A / MS4A6A', 
         GENE == 'EED' ~ 'PICALM / EED', 
         GENE == 'EED' ~ 'PICALM / EED', 
         GENE == "SPDYE3" ~ "PILRA",
         intext == TRUE | microglia_gene == TRUE ~ GENE, 
         ),
       # label = ifelse(!is.na(log_fc), GENE, ""),
       label = ifelse(is.na(label), "", label), 
       label = str_replace_all(label, "e", "ε"),
       microglia_expression_percent = ifelse(is.na(microglia_expression_percent), 0, microglia_expression_percent), 
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
  select(SNP, CHR, BP, GENE, gnomad_maf, OR, architecture, label, microglia_gene, microglia_expression_percent, intext,
         annotation, annotation_impact, dir, effect, size, cytoband, locus, locus_ld, study)
  

## Wrangle Power Curve datasets 
power.dat <- adgwas_power %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or, 
         fill = "fill") %>% 
  filter(study %in% c("Bellenguez")) %>%
  group_by(study) %>%
  fill(or, .direction = "up") %>%
  distinct(maf, or, .keep_all = T) %>%
  ungroup()

## Number of microglia loci
filter(dat.p, microglia_gene == T) %>% group_by(GENE)
filter(dat.p, !str_detect(GENE, "APOE")) %>% group_by(GENE)

# Plotting 
theme.size = 8
geom.text.size = (theme.size - 2) * 0.36
my_colors <- c("grey75", "#377EB8")
names(my_colors) <- c(FALSE, TRUE)

## Sample size plot 
ss.p <- ggplot() + 
  geom_rect(data = tibble(ymax = 870000, ymin = 340000, xmax = 77, xmin = 73), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            color = "grey80", size = 0.25, fill = "grey95"
  ) + 
  geom_point(data = adgwas, aes(x = n_loci, y = size, color = model)) + 
  geom_text(aes(x = 70, y = 530000, label = "Bellenguez"), size = geom.text.size, angle = 90) + 
  geom_text(aes(x = 38, y = 332376, label = "Wightman"), size = geom.text.size, angle = 90, hjust = 0, nudge_y = 0.07) + 
  geom_text(aes(x = 29, y = 455258, label = "Jansen"), size = geom.text.size, angle = 90, hjust = 0, nudge_y = 0.07) +
  geom_text(aes(x = 24, y = 94437, label = "Kunkle"), size = geom.text.size, angle = 90, hjust = 0, nudge_y = 0.07) + 
  geom_text(aes(x = 19, y = 74046, label = "Lambert"), size = geom.text.size, angle = 90, hjust = 0, nudge_y = 0.07) + 
  # geom_text_repel(
  #   # data = adgwas %>% group_by(study) %>% slice(which.max(size)) %>% ungroup(),
  #   data = adgwas %>% mutate(study = as.character(study), study = ifelse(model == 'neff', "", study)),
  #   aes(x = n_loci, y = size, label = study),
  #   vjust = 0,
  #   hjust = 0,
  #   color = 'black',
  #   angle = 90,
  #   # nudge_y = 20000,
  #   size = geom.text.size, 
  #   label.size = NA
# ) +
scale_x_continuous(breaks = c(19, 29, 38, 75)) +
  scale_color_manual(values = c("grey75", "orange"), labels = c("Total", "Effective")) + 
  scale_y_continuous(trans = 'log10', 
                     # labels = scales::comma, 
                     labels = c("50K", "100K", "300K", "1M"), 
                     breaks = c(50000, 100000, 300000, 1000000)) + 
  labs(x = "Number of loci", y = "Sample Size") + 
  # guides(color=guide_legend(override.aes=list(fill=NA))) + 
  theme_light() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        text = element_text(size = theme.size), 
        legend.title = element_blank(), 
        legend.position = c(0.75,0.1), 
        # legend.position = "none", 
        plot.background = element_rect(fill = "transparent", colour = NA_character_), 
        panel.background = element_rect(fill = "transparent", colour = NA_character_), 
        legend.background =  element_rect(colour = "transparent", fill = "transparent"), 
        # legend.key =  element_rect(colour = "transparent", fill = "transparent")
        legend.key = element_blank(), 
        legend.spacing.y = unit(0.0, 'cm'), 
        legend.key.size = unit(0.3, "cm")
        
  ) 

ss.p

## Genetic Architecture plot 
adgwas.p <- ggplot() + 
  geom_rect(data = tibble(ymax = 1.27, ymin = 0.75, xmax = 0.53, xmin = 0.049), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            color = "grey80", size = 0.25, fill = "grey95"
  ) + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_line(data = power.dat, aes(x = maf, y = or), lwd=0.25, color ="orange") + 
  geom_line(data = power.dat, aes(x = maf, y = inv_or), lwd=0.25, color ="orange") + 
  geom_text_repel(
    data = filter(dat.p, dir == "risk" & architecture %nin% c("Common, Low", "Common, Moderate")),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = 0.2, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(dat.p, dir == "protective" & architecture %nin% c("Common, Low", "Common, Moderate")),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = -0.2, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  # geom_point(data = dat.p, 
  #            aes(x = gnomad_maf, 
  #                y = OR, 
  #                fill = microglia_gene, 
  #                color = intext
  #                # size = microglia_expression_percent
  #            ), shape = 21) +
  geom_point(data = dat.p, 
             aes(x = gnomad_maf, 
                 y = OR, 
                 color = microglia_gene, 
                 ), shape = 20) +
  scale_color_manual(values = my_colors,
                     breaks=c("TRUE"),
                     labels=c("Human Microglia\nSignature Gene")
  ) +
  theme_classic() + 
  scale_x_continuous(
    trans='log10',
    breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
    labels = c("5e-1", "0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"),
    # limits = c(0.0001, 0.5)
    ) + 
  scale_y_continuous(
    trans='log',
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25", "0.5", "1", "2", "4", "8"),
    # limits = c(0.2, 13)
    ) + 
  labs(x = "Global Population Minor Alelle Frequency", 
       y = "Odds Ratio - Minor Allele") + 
  coord_cartesian(xlim=c(0.00002, 0.5), ylim = c(0.2, 13)) + 
  theme(
    text = element_text(size = theme.size), 
    legend.title = element_blank(),
    legend.position=c(0.9, 0.95), 
    legend.background = element_rect(linetype="solid"),
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
  )

adgwas.p

## Facet zoom common loci
cutout.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_smooth(data = power.dat, aes(x = maf, y = or), lwd=0.25, color ="orange", se = F) + 
  geom_smooth(data = power.dat, aes(x = maf, y = inv_or), lwd=0.25, color ="orange", se = F) + 
  geom_text_repel(
    data = filter(dat.p, dir == "risk" & architecture %in% c("Common, Low", "Common, Moderate")),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = 0.05, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(dat.p, dir == "protective" & architecture %in% c("Common, Low", "Common, Moderate")),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = -0.05, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_point(data = dat.p, 
             aes(x = gnomad_maf, 
                 y = OR, 
                 color = microglia_gene, 
             ), shape = 20) +
  scale_color_manual(values = my_colors,
                     breaks=c("TRUE"),
                     labels=c("Human Microglia\nSignature Gene")
  ) +
  theme_light() + 
  scale_x_continuous(
    trans='log10',
    breaks = c(0.00001, 0.0001, 0.05, 0.01, 0.1, 0.25, 0.5),
    labels = c("5e-1", "0.0001", "0.05", "0.01", "0.1", "0.25", "0.5"),
    # limits = c(0.0001, 0.5)
  ) + 
  scale_y_continuous(
    trans='log',
    breaks = c(0.8, 0.9, 1, 1.1, 1.2),
    labels = c("0.8", '0.9', "1", "1.1", "1.2"),
    # limits = c(0.2, 13)
  ) + 
  coord_cartesian(xlim=c(0.051, 0.5), ylim = c(0.75, 1.3)) + 
  theme(
    text = element_text(size = theme.size), 
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.background = element_rect(linetype="solid"),
    legend.position = 'none', 
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
  )

cutout.p

## Zoom-in polygon 
positions = tribble(
  ~x, ~y, 
  0.1, 0,
  7.3, 0.98,
  # 6.75, 1.85,
  # 8.65, 1.85, 
  8.72, 0.98, 
  9, 0
)

ribbon.p <- ggplot() +
  geom_segment(aes(x = 1.52, y = 1.84, xend = 2.16, yend = 2.02), 
               size = 0.25, color = 'grey90') + 
  geom_segment(aes(x = 1.52, y = 1.42, xend = 2.15, yend = 0.45), 
               size = 0.25, color = 'grey90') + 
  geom_polygon(data = positions, aes(x = x, y = y), fill = 'grey95') + 
  theme_void() +
  coord_cartesian(xlim=c(0, 9), ylim = c(0, 2)) + 
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA_character_), 
    panel.background = element_rect(fill = "transparent", colour = NA_character_)
  )

ribbon.p

## The empty void 
void.p <- ggplot() + theme_void() + 
  theme(
    # panel.background = element_rect(fill='white'),
    plot.background = element_rect(fill='white', color=NA),
  )

## Combined plot
outfile_neuron.path = '~/Dropbox/Research/PostDoc-MSSM/ADGenetics/results/plots/Neuron_AD_GeneticArchitecture.png'
# tiff(outfile_neuron.path, width = 9, height = 4.5, units = "in", res = 300)
png(outfile_neuron.path, width = 9, height = 4.5, units = "in", res = 600)
pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = ribbon.p,
  x = 0, y = 0,
  width = 9.35, height = 2.18, just = c("left", "top")
)

plotGG(
  plot = void.p,
  x = 2.05, y = 0.5,
  width = 0.15, height = 0.75, just = c("left", "top")
)


plotGG(
  plot = ss.p,
  x = 0, y = 0,
  width = 2, height = 2, just = c("left", "top")
)


plotGG(
  plot = adgwas.p,
  x = 2, y = 0,
  width = 7, height = 2, just = c("left", "top")
)

plotGG(
  plot = cutout.p,
  x = 0.23, y = 2,
  width = 8.77, height = 2.5, just = c("left", "top")
)

pageGuideHide()
dev.off()


## Scater pie plot for experssion levels
cell_exp <- c("glutamatergic_neurons_expression_percent", "gab_aergic_neurons_expression_percent", 
              "astrocytes_expression_percent", "microglia_expression_percent", 
              "oligodendrocytes_expression_percent", "endothelial_cells_expression_percent") 

exp.p <- filter(dat.p, !is.na(microglia_expression_percent))

adexp.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
   geom_text_repel(
    data = filter(exp.p, dir == "risk"),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = 0.15, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(exp.p, dir == "protective"),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = -0.25, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  scatterpie::geom_scatterpie(aes(x=gnomad_maf, y=OR), cols=cell_exp, color=NA, data=exp.p) + 
  theme_classic() + 
  coord_equal() +
  labs(x = "Global Population Minor Alelle Frequency", 
       y = "Odds Ratio - Minor Allele") + 
  scale_fill_manual(name = "Single Cell Experession", 
                    labels = c("Glutamatergic Neurons", "Gabaergic Neurons", "Astrocytes", 
                               "Microglia", "Oligodendrocytes", "Endothelial cells"),
                    values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF")) +
  # scale_colour_manual(values = use_col) + 
  guides(color=guide_legend(ncol=1)) + 
  theme(
    text = element_text(size = theme.size), 
    legend.title = element_blank(),
    legend.position= "bottom", 
    legend.background = element_rect(linetype="solid")
    # legend.position = 'none'
  )

ggsave("~/Downloads/AD_expression.png", plot = adexp.p, 
       width = 4.5, height = 4.5, units = "in")


















































