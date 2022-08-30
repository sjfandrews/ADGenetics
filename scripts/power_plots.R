library(tidyverse)

ad_loci <- read_csv("/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/Neurogenomics/intermediate/ad_loci.csv")

adgwas <- tribble(
  ~study, ~year, ~neff, ~n, ~n_pca, ~n_ca, ~n_cn, ~n_loci, ~ancestry, ~cohorts, ~notes,
  "Lambert", 2013, NA, 74046, NA, 17008, 37154, 19, NA, NA, NA,
  "Kunkle", 2019, NA, 94437, NA, 35274, 59163, 24, NA, NA, NA, 
  "Marioni", 2018, NA, 377012, NA, (27696 + 14338 + 25580), (37154 + 245941), 26, NA, NA, NA,
  "Jansen", 2019, NA, 455258, NA, (24087 + 47793), (55058 + 328320), 29, NA, NA, NA,
  "Wightman", 2021, NA, 1126563, NA, 90338, 1036225, 38, NA, NA, NA,
  "Bellenguez", 2022, NA, 528781, NA, (39106 + 46828 + 25392), (401577 + 276086), 75, NA, NA, NA,
  # "Lake", 2023, NA, 644188, NA, (54233 + 46828), 543127, NA, NA, NA, NA,
  "Future1", 2050, NA, 500000, NA, 250000, 250000, NA, NA, NA, NA,
  "Future2", 2100, NA, 1000000, NA, 500000, 500000, NA, NA, NA, NA,
) %>% 
  mutate(
  	effN = 4 / ((1/n_ca) + (1/n_cn)), 
	phi = n_ca/n,
	neff = n*phi*(1-phi)
  )

################# Draw curves indicating the effect size needed across different MAF threshols

# power wanted
pw.thresh = 0.8 
# significance threshold (aka alpha)
p.threshold = 5e-8
# calculate the chi-square value corresponding to significance threshold defined in p.threshold
q = qchisq(p.threshold, df = 1, lower = F) 

# Sequence of frequencies from min MAF to 0.5
f = c(seq(0.0000001, 0.01, length = 1000), seq(0.01, 0.5, length = 500)) 
# Sequence of effect sizes from min beta to max beta
b = seq(0, 4, length = 1500)    


out <- list(NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL) %>% 
  magrittr::set_names(adgwas$study)

for(z in 1:nrow(adgwas)){
  STUDY = adgwas$study[z]
  eff.N = adgwas %>% filter(study == STUDY) %>% pull(neff)
  b.for.f = rep(NA, length(b))
  
  for(i in 1:length(b)){ 
    if(i%%100 == 0){message(STUDY, ": ", i)}
    # Calculate power at this allele frequency, across a range of effect sizes (b)
    pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
    # Calculate what is the minimum b needed to read pw.thres
    b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
  }
  
  out[[z]] = tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
    arrange(maf) %>%
    mutate(inv_or = 1/or,
           study = STUDY)
  
}

power_gwas <- bind_rows(out) %>%
  mutate(
    study = fct_relevel(study, "Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez", "Future1", "Future2")
  )

write_csv(power_gwas, "/Users/sheaandrews/Dropbox/Research/PostDoc-MSSM/ADGenetics/intermediate/ad_power_matti.csv")


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



# --------------- For BINARY: numbers of cases (col1) and controls (col2) to calculate effective sample size
eff.N = adgwas %>% filter(study == "Bellenguez") %>% pull(neff)
# create null variable b.for.f
b.for.f = rep(NA, length(b))

for(i in 1:length(b)){ 
  # Calculate power at this allele frequency, across a range of effect sizes (b)
  pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
  # Calculate what is the minimum b needed to read pw.thres
  b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
}

power_bellenguez <- tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or,
         study = "Bellenguez")

#### Wightman
eff.N = adgwas %>% filter(study == "Wightman") %>% pull(neff)
b.for.f = rep(NA, length(b))

for(i in 1:length(b)){ 
  pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
  b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
}

power_wightman <- tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or,
         study = "Wightman")

#### Jansen
eff.N = adgwas %>% filter(study == "Jansen") %>% pull(neff)
b.for.f = rep(NA, length(b))

for(i in 1:length(b)){ 
  pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
  b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
}

power_jansen <- tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or,
         study = "Jansen")

#### Marioni
eff.N = adgwas %>% filter(study == "Marioni") %>% pull(neff)
b.for.f = rep(NA, length(b))

for(i in 1:length(b)){ 
  pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
  b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
}

power_marioni <- tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or,
         study = "Marioni")

#### Kunkle
eff.N = adgwas %>% filter(study == "Kunkle") %>% pull(neff)
b.for.f = rep(NA, length(b))

for(i in 1:length(b)){ 
  pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
  b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
}

power_kunkle <- tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or,
         study = "Kunkle")

#### Lambert
eff.N = adgwas %>% filter(study == "Lambert") %>% pull(neff)
b.for.f = rep(NA, length(b))

for(i in 1:length(b)){ 
  pwr_bin = pchisq(q, df = 1, ncp=2*f[i]*(1-f[i])*b^2*eff.N, lower=F)
  b.for.f[i] = b[min( which(pwr_bin > pw.thresh))]
}

power_lambert <- tibble(maf = f, obs.b = b, beta = b.for.f, or = exp(b.for.f)) %>%
  arrange(maf) %>%
  mutate(inv_or = 1/or,
         study = "Lambert")

#### 
power_gwas <- bind_rows(
  power_lambert, 
  power_kunkle, 
  power_marioni, 
  power_jansen, 
  power_wightman, 
  power_bellenguez
) %>%
  mutate(
    study = fct_relevel(study, "Lambert", "Marioni", "Jansen", "Kunkle", "Wightman", "Bellenguez")
  )
















