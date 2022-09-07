dat_loci_power <- ad_loci %>% 
  left_join(adgwas, by = "study") %>%
  filter(!str_detect(SNP, "APOE")) %>%
  filter(study != "Jonsson") %>%
  mutate( # TOFIX 
    OR = ifelse(SNP == "rs139643391", exp(-0.0619), OR),
    gnomad_maf = ifelse(SNP == "rs149080927", 0.3829, gnomad_maf),
    dir = ifelse(OR > 1, "risk", "protective"), 
    effect = ifelse(OR > 1, OR, 1/OR), 
    study = fct_relevel(study, "Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez"), 
  ) 

ggplot(power_gwas, aes(y = or, x = maf, color = study)) + 
  geom_line(lwd=0.25) + 
  scale_x_continuous(trans='log10', breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.00001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  scale_color_manual(values = c("Lambert" = "#1F78B4",
                                "Kunkle" = "#33A02C",
                                "Jansen" = "#E31A1C",
                                "Marioni" = "#6A3D9A", 
                                "Wightman" = "#FF7F00", 
                                "Bellenguez" = "#A65628")) + 
  theme_light()

ggplot(dat_loci_power, aes(x = gnomad_maf, y = effect, color = study)) + 
  geom_line(data = power_gwas, aes(y = or, x = maf, color = study), lwd=0.25) + 
  geom_point(alpha = 0.5) + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  scale_color_manual(values = c("Lambert" = "#1F78B4",
                                "Kunkle" = "#33A02C",
                                "Jansen" = "#E31A1C",
                                "Marioni" = "#6A3D9A", 
                                "Wightman" = "#FF7F00", 
                                "Bellenguez" = "#A65628")) +
  facet_wrap(vars(study), nrow = 2) + 
  theme_light()

gwas_power.p <- ggplot(dat_loci_power, aes(x = gnomad_maf, y = effect, color = study)) + 
  geom_line(data = power_gwas, aes(y = or, x = maf, color = study), lwd=0.25) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.5),
                     labels = c("1e-4", "0.001", "0.01", "0.1",  "0.5"), 
                     limits = c(0.0001, 0.5)) + 
  scale_y_continuous(trans='log10') + 
  scale_color_manual(values = c("Lambert" = "#1F78B4",
                                "Kunkle" = "#33A02C",
                                "Jansen" = "#E31A1C",
                                "Marioni" = "#6A3D9A", 
                                "Wightman" = "#FF7F00", 
                                "Bellenguez" = "#A65628")) +
  facet_wrap(vars(study), nrow = 2) +
  labs(x = "Global Population Minor Alelle Frequency", y = "Odds Ratio - Minor Allele") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        text = element_text(size = 8), 
        legend.position = "none", 
        strip.background = element_blank())
gwas_power.p

ggsave("~/Downloads/AD_gwas_power.png", plot = gwas_power.p + theme(legend.position = "none"), 
       width = 6, height = 3, units = "in")


# Sample size vs no. loci Plots
adgwas.p <- adgwas %>%
  mutate(label = paste(study, year, sep = ", ")) %>%
  select(label, n, effN, n_loci) %>%
  slice(1:6) %>%
  pivot_longer(c("n", "effN"), names_to = "model", values_to = "sample_size") %>%
  mutate(
  	model = fct_relevel(model, "n", "effN"),
  	model = fct_recode(model, "Total n" = "n", "Effective n" = "effN")
  )

ad_nloci.p <- ggplot(adgwas.p, aes(x = n_loci, y = sample_size, label = label, color = label)) + 
  geom_point()  + 
  scale_color_manual(values = c("Lambert, 2013" = "#1F78B4",
                                "Kunkle, 2019" = "#33A02C",
                                "Jansen, 2019" = "#E31A1C",
                                "Marioni, 2018" = "#6A3D9A", 
                                "Wightman, 2021" = "#FF7F00", 
                                "Bellenguez, 2022" = "#A65628")) + 
  ggrepel::geom_text_repel(size = 2, color = 'black') + 
  facet_wrap(vars(model), nrow = 2, scales = "free_y") + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y = "Sample Size", x = "Loci discovered") + 
  theme_bw() + 
  theme(text = element_text(size = 8), 
        aspect.ratio = 1,
        strip.background = element_blank(),
        legend.position = "none", 
        axis.text.y = element_text(angle = 90, hjust = 0.5) 
  )

ad_nloci.p

## Plot Gardner
png("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/plots/AD_gwas_power.png", 
    width = 9, height = 4.5, units = "in", res = 300)
pageCreate(width = 9, height = 4.5, default.units = "inches")

plotGG(
  plot = ad_nloci.p ,
  x = 0, y = 0,
  width = 2.25, height = 4.5, just = c("left", "top")
)

plotGG(
  plot = gwas_power.p,
  x = 2, y = 0,
  width = 6.75, height = 4.5, just = c("left", "top")
)

pageGuideHide()
dev.off()













## Effective sample size
p.neff <- ggplot(adgwas.p, aes(x = n_loci, y = neff, label = label, color = label)) + 
  geom_point()  + 
  scale_color_manual(values = c("Lambert, 2013" = "#1F78B4",
                                "Kunkle, 2019" = "#33A02C",
                                "Jansen, 2019" = "#E31A1C",
                                "Marioni, 2018" = "#6A3D9A", 
                                "Wightman, 2021" = "#FF7F00", 
                                "Bellenguez, 2022" = "#A65628")) + 
  ggrepel::geom_label_repel(size = 2, color = 'black') + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y = "Effective Sample Size", x = "Loci discovered") + 
  theme_bw() + 
  theme(text = element_text(size = 8), 
        aspect.ratio = 1,
        legend.position = "none")

## Total sample size
p.n <- ggplot(adgwas.p, aes(x = n_loci, y = n, label = label, color = label)) + 
  geom_point() + 
  scale_color_manual(values = c("Lambert, 2013" = "#1F78B4",
                                "Kunkle, 2019" = "#33A02C",
                                "Jansen, 2019" = "#E31A1C",
                                "Marioni, 2018" = "#6A3D9A", 
                                "Wightman, 2021" = "#FF7F00", 
                                "Bellenguez, 2022" = "#A65628")) + 
  ggrepel::geom_label_repel(size = 2, color = 'black') + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y = "Total Sample Size", x = "Loci discovered") + 
  theme_bw() + 
  theme(text = element_text(size = 8), 
        aspect.ratio = 1,
        legend.position = "none")

p.loci_n <- cowplot::plot_grid(
  p.n, p.neff, ncol = 1
)

ggsave(plot = p.loci_n, "~/Dropbox/Research/PostDoc-MSSM/Neurogenomics/plots/loci_samplesize.png", 
       units = "cm", height = 13, width = 6
)


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

### Plot loci by year discovered with corresponding poewr lines
adgwas_variants <- adgwas_loci %>% 
  filter(study %nin% c("Jonsson", "Reiman")) %>%
  arrange(locus_ld) %>% 
  mutate(locus_ld = ifelse(is.na(locus_ld), row(.), locus_ld), 
         effect = ifelse(OR > 1, OR, 1/OR), 
         study = fct_relevel(study, "Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez")
         ) %>%
  left_join(adgwas_meta, by = "study") %>%
  group_by(locus_ld) %>%
  arrange(year) %>%
  slice(1) %>%
  ungroup() %>% 
  select(SNP, CHR, BP, GENE, gnomad_maf, OR, annotation, annotation_impact, effect, cytoband, locus, locus_ld, study, year) 

power.dat <- adgwas_power %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or, 
         fill = "fill", 
         study = fct_relevel(study, "Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez", "Future2")) %>% 
  filter(study %nin% c("Future1")) %>%
  # filter(study %in% c("Bellenguez", "Future1", "Future2")) %>%
  group_by(study) %>%
  fill(or, .direction = "up") %>%
  distinct(maf, or, .keep_all = T) %>%
  ungroup()

test <- ggplot(power.dat, aes(x = maf, y = or, color = study), lwd=0.25) + 
  geom_line(alpha = 0.5) + 
  geom_point(data = adgwas_variants, aes(x = gnomad_maf, y = effect, color = study), size = 0.5) + 
  scale_color_brewer(palette = "Set2") + 
  scale_x_continuous(trans='log10',
                     breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("5e-1", "0.0001", "0.001", "0.01", "0.1", "0.25", "0.5")) + 
  scale_y_continuous(trans='log', 
                     breaks=c(0, 1, 2, 4, 8, 12, 16, 20), 
                     labels = c("0", "1", "2", "4", "8", "12", "16", "20")) +
  coord_cartesian(xlim=c(0.00002, 0.5), ylim = c(0.9, 20)) + 
  theme_classic()

ggsave(plot = test, "sandbox/plots/loci_year_discoverd.png", 
       units = "in", height = 6, width = 6
)

