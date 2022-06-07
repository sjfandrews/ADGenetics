library(tidyverse)
library(glue)
library(ggrepel)
library(plotgardener)
library(GGally)
`%nin%` = negate(`%in%`)

setwd("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/Neurogenomics/")

ad_loci <- read_csv("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/Neurogenomics/intermediate/ad_loci.csv")

ad_genes <- select(ad_loci, study, locus, CHR, BP, P, gencode_gene, gvc_gene) %>%
  mutate(gene = ifelse(is.na(gvc_gene), gencode_gene, gvc_gene), 
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
  
upset.dat <- ad_loci %>%
  filter(study %nin% c("Jonsson", "Reiman")) %>%
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
         y_text = 1
         )

trem2 <- filter(upset.dat, CHR ==  6, BP == 41129207)  

upset.sig <- upset.dat %>%
  group_by(study, locus) %>%
  summarize(n_sig = sum(sig), gene = first(gene)) %>%
  ungroup()

upset_lines <- upset.dat %>%
  filter(sig == T) %>%
  group_by(locus) %>%
  slice(which.min(study), which.max(study)) %>%
  select(study, locus, gene, label) %>%
  mutate(pos = c("start", "end")) %>%
  ungroup() %>%
  pivot_wider(names_from = pos, values_from = study)

## Plots 
theme.size = 8
geom.text.size = (theme.size - 2) * 0.36

upset.p <- ggplot(upset.dat, aes(y = study, x = gene, color = sig)) + 
  geom_stripped_rows(colour = "NA") +
  # geom_stripped_cols(colour = "NA") +
  geom_point() +
  # geom_point(data = trem2, aes(y = study, x = gene), color = "red") +
  geom_segment(data = upset_lines, aes(x = gene, y = start, xend = gene, yend = end, colour = "segment")) +
  theme_light() + 
  scale_color_manual(values = c("grey75", "black", "black")) + 
  labs(x = "locus") + 
  theme(
    text = element_text(size = theme.size, family = "sans"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0), 
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

upset.p

upset.p <- ggplot(upset.sig, aes(y = study, x = gene, color = as.factor(n_sig))) + 
  geom_stripped_rows(colour = "NA") +
  # geom_stripped_cols(colour = "NA") +
  geom_point(data = filter(upset.sig, n_sig == 0), aes(y = study, x = gene, color = as.factor(n_sig))) +
  geom_segment(data = upset_lines, aes(x = gene, y = start, xend = gene, yend = end, colour = "segment")) +
  geom_point(data = filter(upset.sig, n_sig > 0), aes(y = study, x = gene, color = as.factor(n_sig))) +
  # geom_point(data = trem2, aes(y = study, x = gene), color = "red") +
  scale_color_manual(values = c("grey75", "black", "#fcae91", "#fb6a4a", "black", "#cb181d")) + 
  theme_light() + 
  labs(x = "locus") + 
  theme(
    text = element_text(size = theme.size, family = "sans"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

upset.p


locus.p <- ggplot(upset.dat, aes(y = y_text, x = label)) + 
  geom_text(data = upset.dat, aes(y = y_text, x = label, label = locus), 
            angle = 90, vjust = 0.5, hjust = 1, size = geom.text.size, 
            family = "sans", color = "grey30") +
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

test_rect <- tribble(
  ~study, ~xmin, ~xmax, ~ymin, ~ymax, 
  "Lambert", -Inf, Inf, -Inf, 1.5,
  "Jansen",  -Inf, Inf, 2.5, 3.5,
  "Wightman", -Inf, Inf, 4.5, 5.5,
)

test %>% 
  group_by(study) %>%
  summarise(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y))

ggplot(aes(x = c(1:86))) + 
  geom_rect(data = test_rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill = "grey75", color = "grey75", alpha=0.5, stat="identity") +
  geom_point() +
  geom_segment(data = test_lines, aes(x = label, y = start, xend = label, yend = end, colour = "segment"))
  

png("~/Dropbox/Research/PostDoc-MSSM/Neurogenomics/plots/AD_upset2.png", width = 9, height = 4.5, units = "in", res = 300)
pageCreate(width = 9, height = 3.5, default.units = "inches")

plotGG(
  plot = locus.p,
  x = 0.46, y = 1.55,
  width = 8.54, height = 2.5, just = c("left", "top")
)

plotGG(
  plot = upset.p,
  x = 0, y = 0,
  width = 9, height = 2, just = c("left", "top")
)

pageGuideHide()
dev.off()





















  