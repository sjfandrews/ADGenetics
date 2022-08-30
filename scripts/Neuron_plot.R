
## From Bellenguez
microglia_gene_list <- c("SPI1", "MAF", "MEF2C", "NR1H2", "NR1H3", "TREM2", "BLNK", "MS4A4A", "MS4A6A", 
                         "CD33", "SIGLEC11", "INPP5D/SHIP1", "NCK2", "ABI3", "PICALM", "EED", "CD2AP", 
                         "BIN1", "SORL1", "RAB10", "SPPL2A", "PGRN", "ZYX", "RIN3", "AP4E1", "AP4M1", 
                         "TP3INP1", "APOE", "CLU", "LACTB", "OPA3")

st20 <- readxl::read_xlsx("resources/Bellenguez/41588_2022_1024_MOESM4_ESM.xlsx", sheet = 21, skip = 2) %>% 
  janitor::clean_names() %>%
  select(gene, locus, ends_with("percent")) %>%
  rename(locus_bellenguez = locus) %>%
  mutate_at(vars(ends_with("percent")), as.numeric)

## Gene label appears in text
intext_list <- c("IDUA", "JAZF1", "TSPAN14", "BLNK", "PLEKHA1", "CTSH", "SIGLEC11", 
            "CR1", "BIN1", "INPP5D", "HLA-DQA1", "TREM2", "EPDR1", "PTK2B", 
            "USP6NL", "EED", "PLCG2", "ABI3", "TSPOAP1", "SPPL2A", "SPI1", "ABCA1", 
            "MS4A4A", "SORL1", "MAF"
)

gosselin.raw <- readxl::read_xlsx("resources/gosselin/aal3222_gosselin_tables2.xlsx", sheet = 2) %>% 
  janitor::clean_names()

gosselin <- gosselin.raw %>% 
  mutate(microglia_gene = TRUE)


dat.p <- dat_loci %>% 
  filter(study %in% c("Bellenguez", "Reiman"))  %>%
  left_join(gosselin, by = c("GENE" = "gene_name")) %>%
  left_join(st20, by = c("GENE" = "gene")) %>%
  mutate(size = ifelse(effect < 1.2, 1, effect), 
         microglia_gene = ifelse(is.na(microglia_gene), FALSE, microglia_gene),
         intext = ifelse(GENE %in% intext_list | str_detect(GENE, "APOE"), TRUE, FALSE), 
         # microglia_gene = case_when(
         #   GENE %in% microglia_gene_list ~ TRUE,
         #   TRUE ~ FALSE),
         label = case_when(
           annotation_impact %in% c("HIGH", "LOW", "MODERATE") & gnomad_maf < 0.05  & microglia_gene == TRUE ~ glue("{gene_name} {hgvs_p_new}"),
           str_detect(GENE, "APOE") ~ GENE,
           intext == TRUE | microglia_gene == TRUE ~ GENE),
         # label = ifelse(!is.na(log_fc), GENE, ""),
         label = ifelse(is.na(label), "", label), 
         label = str_replace_all(label, "ε", "e"),
         microglia_expression_percent = ifelse(is.na(microglia_expression_percent), 0, microglia_expression_percent)
  ) %>% 
  select(SNP, CHR, BP, GENE, gnomad_maf, OR, label, architecture, architecture2, microglia_gene, microglia_expression_percent, intext,
         annotation, annotation_impact, dir, effect, size, cytoband, locus, locus_ld, study)
  
 
## Number of microglia loci
filter(dat.p, microglia_gene == T) %>% group_by(GENE)
filter(dat.p, !str_detect(GENE, "APOE")) %>% group_by(GENE)

theme.size = 8
geom.text.size = (theme.size - 2) * 0.36
my_colors <- c("grey75", "#377EB8")
names(my_colors) <- c(FALSE, TRUE)


adgwas.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_text_repel(
    data = filter(dat.p, dir == "risk"),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = 0.15, seed = 333,
    segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
    color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(
    data = filter(dat.p, dir == "protective"),
    aes(x = gnomad_maf, y = OR, label = label),
    max.overlaps = Inf, nudge_y = -0.15, seed = 333,
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
    # trans='log10',
    breaks = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
    limits = c(0.001, 0.5)) + 
  scale_y_continuous(
    trans='log',
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25", "0.5", "1", "2", "4", "8"),
    limits = c(0.2, 13)) + 
  labs(x = "Global Population Minor Alelle Frequency", 
       y = "Odds Ratio - Minor Allele") + 
  theme(
    text = element_text(size = theme.size), 
    legend.title = element_blank(),
    legend.position=c(0.9, 0.95), 
    legend.background = element_rect(linetype="solid")
    # legend.position = 'none'
  )

adgwas.p

ggsave("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADgenetics/plots/AD_microglia_gosslin.png", plot = adgwas.p, 
       width = 7.4, height = 4.5, units = "in")
ggsave("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADgenetics/plots/AD_microglia_gosslin.pdf", plot = adgwas.p, 
       width = 7.4, height = 4.5, units = "in")


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


















































