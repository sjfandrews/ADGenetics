# Estimate Power curves 
## https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS3.html
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

## Snakemake 
### Input
adgwas_meta.path =  "intermediate/gwas_metadata.csv"

adgwas_meta.path = snakemake@input[['meta']]

### Output
outfile = snakemake@output[['power_curves']]

## Import data 
adgwas_meta.raw <- read_csv(adgwas_meta.path)

adgwas_meta <- adgwas_meta.raw %>%
  filter(study %in% c("Lambert", "Kunkle", "Marioni", "Jansen", "Wightman", "Bellenguez")) %>%
  add_row(study = "Future1", year = 2027, n = 500000, n_ca = 250000, n_cn = 250000) %>%
  add_row(study = "Future2", year = 2032, n = 1000000, n_ca = 500000, n_cn = 500000)  %>% 
  mutate(
    effN = 4 / ((1/n_ca) + (1/n_cn)), 
    phi = n_ca/n,
    neff = n*phi*(1-phi)
  )

### Draw curves indicating the effect size needed across different MAF thresholds
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


out <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) %>% 
  magrittr::set_names(adgwas_meta$study)

message("\nCalculating Power....\n")

for(z in 1:nrow(adgwas_meta)){
  STUDY = adgwas_meta$study[z]
  eff.N = adgwas_meta %>% filter(study == STUDY) %>% pull(neff)
  b.for.f = rep(NA, length(b))
  
  for(i in 1:length(b)){ 
    if(i%%100 == 0){message("\t", STUDY, ": ", i)}
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

write_csv(power_gwas, outfile)


















