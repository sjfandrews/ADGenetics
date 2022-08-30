# Upset plot of AD Loci 
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
library(GGally)
`%nin%` = negate(`%in%`)

## Snakemake 
### Input
ad_loci.path = "results/adgwas_loci.csv"
adgwas_meta.path =  "intermediate/gwas_metadata.csv"

ad_loci.path = snakemake@input[['loci']]
adgwas_meta.path = snakemake@input[['meta']]

### Output
upset.path = snakemake@output[['upset']]

## Import Datasets
ad_loci <- read_csv(ad_loci.path)
adgwas.raw <- read_csv(adgwas_meta.path)

## Wrangle data 
### AD GWAS meta data 
adgwas <- adgwas.raw %>% 
  filter(study %in% c("Lambert", "Kunkle", "Marioni", "Jansen", "Wightman", "Bellenguez")) %>%
  mutate(
    study = fct_relevel(study, c("Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez")), 
    n = n - neff,
    ) %>%
  select(study, year, neff, n) %>%
  pivot_longer(c("neff", "n"), names_to = "model", values_to = "size") 

### AD Loci + Gene lables
ad_genes <- select(ad_loci, study, locus, CHR, BP, P, GENE, gencode_gene, gvc_gene) %>%
  filter(study %in% c("Lambert", "Kunkle", "Marioni", "Jansen", "Wightman", "Bellenguez")) %>%
  mutate(gene = GENE,
         gene = ifelse(is.na(gvc_gene), gencode_gene, gvc_gene), 
         y = case_when(
           study == "Lambert" ~ 1, 
           study == "Marioni" ~ 2, 
           study == "Jansen" ~ 3, 
           study == "Kunkle" ~ 4, 
           study == "Wightman" ~ 5, 
           study == "Bellenguez" ~ 6, 
         )) %>%
  group_by(locus) %>%
  arrange(-y) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(CHR, BP) %>%
  mutate(x = 1:n()) %>%
  select(locus, gene, x) 

### AD loci with significance   
upset.dat <- ad_loci %>%
  filter(study %in% c("Lambert", "Kunkle", "Marioni", "Jansen", "Wightman", "Bellenguez")) %>%
  filter(!is.na(locus)) %>%
  expand(study, locus) %>%
  left_join(select(ad_loci, study, locus, locus_ld, CHR, BP, P)) %>%
  left_join(ad_genes, by = "locus") %>%
  arrange(CHR, BP) %>%
  mutate(label = glue("{locus} ({gene})"),
         sig = ifelse(is.na(P), FALSE, TRUE), 
         sig = ifelse(gene == "APOE", TRUE, sig), 
         # locus_ld = ifelse(is.na(locus_ld), 1, locus_ld), 
         study = fct_relevel(study, "Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez"), 
         locus = fct_inorder(locus),
         label = fct_inorder(label),
         gene = fct_inorder(gene),
         y = case_when(
           study == "Lambert" ~ 1, 
           study == "Marioni" ~ 2, 
           study == "Jansen" ~ 3, 
           study == "Kunkle" ~ 4, 
           study == "Wightman" ~ 5, 
           study == "Bellenguez" ~ 6, 
         ), 
         y_text = 1,
         )  %>%
  separate(locus, into = c("chr", "start", "end")) %>%
  mutate(start = round(as.numeric(start) / 1000000, 2), 
         end = round(as.numeric(end) / 1000000, 2), 
         locus = glue("{chr}:{start}-{end}")
  )
  
### Significant loci
upset.sig <- upset.dat %>%
  group_by(study, locus) %>%
  summarize(n_sig = sum(sig), gene = first(gene)) %>%
  ungroup() %>%
  # filter(n_sig > 1) %>%
  left_join(upset.dat)

### Joining lines
upset_lines <- upset.dat %>%
  filter(sig == T) %>%
  group_by(locus) %>%
  slice(which.min(study), which.max(study)) %>%
  select(study, locus, x, label) %>%
  mutate(pos = c("start", "end")) %>%
  ungroup() %>%
  pivot_wider(names_from = pos, values_from = study)

## Plots
theme.size = 8
geom.text.size = (theme.size - 2) * 0.36

### Sample size bar graph
adgwas_ss.p <- ggplot(adgwas, aes(x = size, y = study, fill = model)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("black", "orange"), labels = c("Total", "Effective")) +
  labs(x = "Sample Size")  + 
  theme_classic() +  
  scale_x_reverse(breaks = c(300000, 600000, 900000), labels = c("300k", "600k", "900k")) + 
  # guides(fill = guide_legend(override.aes = list(size = 4))) + 
  theme(
    text = element_text(size = theme.size), 
    axis.title.y=element_blank(), 
    axis.text.y=element_blank(), 
    axis.ticks.y=element_blank(), 
    legend.title = element_blank(), 
    legend.position=c(0.2, 0.2), 
    legend.key.size = unit(0.25, 'cm'),
    axis.line.y = element_blank()
  )

adgwas_ss.p

### Upset plot showing common loci
upset.p <- ggplot(upset.dat, aes(y = study, x = x)) + 
  geom_point(data = upset.sig %>% filter(n_sig == 0), aes(y = study, x = x), color = "grey90") + 
  geom_segment(data = upset_lines, aes(x = x, y = start, xend = x, yend = end, colour = "segment")) +
  geom_stripped_rows(colour = "NA") +
  # geom_point(data = upset.sig %>% filter(n_sig > 1), aes(y = study, x = x + 0.2), color = "#ef3b2c") + 
  geom_point(data = upset.sig %>% filter(n_sig > 1), aes(y = study, x = x), 
             shape = 21, fill = "#ef3b2c", color = 'black', stroke = 1, size = 0.9) + 
  geom_point(data = upset.sig %>% filter(n_sig == 1), aes(y = study, x = x), color = "black") +
  # geom_point() +
  theme_light() + 
  scale_color_manual(values = c("black", "black", "black")) + 
  labs(x = "locus") + 
  scale_x_continuous(breaks = unique(upset.dat$x),
                     labels = unique(upset.dat$gene),
                     ) + 
  coord_cartesian(xlim=c(4,78)) +
  theme(
    text = element_text(
      size = theme.size, 
      # family = "sans"
      ),
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0,
    #                            face=ifelse(is.na(ad_genes$gvc_gene),"plain","bold")),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
  )

upset.p

### Locus labels 
locus.p <- ggplot(upset.dat, aes(y = y_text, x = label)) + 
  geom_text(data = upset.dat, aes(y = y_text, x = label, label = locus), 
            angle = 90, vjust = 0.5, hjust = 1, size = geom.text.size, 
            # family = "sans", 
            color = "grey30") +
  # geom_text(data = test, aes(y = y_text+0.1, x = label, label = gene),
  #           angle = 90, vjust = 0.5, hjust = 0, size = geom.text.size,
  #           family = "sans", color = "grey30") +
  # coord_cartesian(ylim=c(0.98,1)) +
  coord_fixed(ratio = 500, ylim=c(0.97,1)) + 
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

locus.p

## Export Plot
png(upset.path, width = 9, height = 4.5, units = "in", res = 300)
pageCreate(width = 9, height = 3.5, default.units = "inches")

plotGG(
  plot = adgwas_ss.p,
  x = 0, y = 0,
  width = 1.6, height = 1.75, just = c("left", "top")
)


plotGG(
  plot = locus.p,
  x = 1.98, y = 1.44,
  width = 7, height = 2.5, just = c("left", "top")
)


plotGG(
  plot = upset.p,
  x = 1.5, y = 0,
  width = 7.5, height = 2, just = c("left", "top")
)

pageGuideHide()
dev.off()


















  