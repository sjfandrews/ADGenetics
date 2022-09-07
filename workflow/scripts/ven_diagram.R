library(nVennR)
library(tidyverse)

# lambert
78 + 353 + 68 + 117 + 66 + 462 + 171 

sample_size <- read_csv("~/Dropbox/Research/PostDoc-MSSM/ADGenetics/resources/ad_sample_sizes.csv")

sample_size %>%
  rowwise() %>%
  mutate(wightman = sum(wightman_cases, wightman_controls, na.rm = T), 
         bellenguez = sum(bellenguez_cases, bellenguez_controls, na.rm = T), 
         wightman_samps = paste0(cohort, seq(1:wightman))) %>%
  ungroup() %>%
  
  group_by(study) %>%
  
  summarise(tot = sum(tot))


test <- sample_size %>%
  filter(., str_detect(cohort, "CHARGE-")) %>%
  rowwise() %>%
  mutate(wightman = sum(wightman_cases, wightman_controls, na.rm = T), 
         bellenguez = sum(bellenguez_cases, bellenguez_controls, na.rm = T), 
         lambert = sum(lambert_cases, lambert_controls, na.rm = T), 
         marioni = sum(marioni_cases, marioni_controls, na.rm = T), 
         wightman_samps = list(paste(cohort, seq(1:wightman), sep = "_")), 
         bellenguez_samps = list(paste(cohort, seq(1:bellenguez), sep = "_")), 
         lambert_samps = list(paste(cohort, seq(1:lambert), sep = "_")),
         marioni_samps = list(paste(cohort, seq(1:marioni), sep = "_"))
         )

sample_size %>%
  filter(., str_detect(cohort, "CHARGE-")) %>%
  summarise(
    lambert_cases = sum(lambert_cases, na.rm = T), 
    lambert_controls = sum(lambert_controls, na.rm = T),
    wightman_cases = sum(wightman_cases, na.rm = T), 
    wightman_controls = sum(wightman_controls, na.rm = T),
    bellenguez_cases = sum(bellenguez_cases, na.rm = T), 
    bellenguez_controls = sum(bellenguez_controls, na.rm = T)
  )


# bell only
bell_prefix <- "bell"
bell_suffix <- seq(1:51968)
bell = paste(bell_prefix, bell_suffix, sep = "")

# wight only
wight_prefix <- "wight"
wight_suffix <- seq(1:611917)
wight = paste(wight_prefix, wight_suffix, sep = "")

# both
both_prefix <- "both"
both_suffix <- seq(1:1251777)
both = paste(both_prefix, both_suffix, sep = "")


Bellenguez = c(bell, both)
Wightman = c(wight, both)

myV <- plotVenn(list(
  'Lambert' = test$lambert_samps,
  'Marioni' = test$marioni_samps,
  'Bellenguez' = test$bellenguez_samps,
  'Wightman' = test$wightman_samps), 
                nCycles = 6000)

showSVG(myV, outFile = '~/Downloads/sampleoverlap.svg', 
        opacity = 0.1, borderWidth = 3, labelRegions = F, fontScale = 2.2, 
        setColors=c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A"))
fixsvg('~/Downloads/sampleoverlap.svg')


#setwd('/Users/shea/Dropbox/Research/PostDoc-MSSM/8_ADGWAS')
## -------------------------------------------------------------------------------- ##
##        Mock Sample IDs for Sample overlap Venn

fixsvg <- function(svgname) {
  decfind <- function(svglines) {
    stringr::str_detect(svglines, "<svg") & !stringr::str_detect(svglines, "www.w3.org/1999/xlink")
  }
  svg <- readLines(svgname)
  declaration <- svg[decfind(svg)]
  strt <- substring(declaration, 1, nchar(declaration) - 1)
  svg[decfind(svg)] <- paste0(
    strt, ' xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">')
  writeLines(svg, svgname)
}

breaksvg <- function(svgname) {
  legfind <- function(svglines) {
    stringr::str_detect(svglines, "class=\"([pq][0-9]|legend)\"") |
      stringr::str_detect(svglines, "www.w3.org/1999/xlink")
  }
  svg <- readLines(svgname)
  declaration <- svg[decfind(svg)]
  strt <- substring(declaration, 1, nchar(declaration) - 1)
  svg[decfind(svg)] <- paste0(
    strt, ' xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">')
  writeLines(svg, svgname)
}

## -------------------------------------------------------------------------------- ##
##        Mock Sample IDs for Sample overlap Venn

# n1234
n1234_prefix <- "abcd"
n1234_suffix <- seq(1:54162)
n1234 = paste(n1234_prefix, n1234_suffix, sep = "")

# n124
n124_prefix <- "abd"
n124_suffix <- seq(1:19884)
n124 = paste(n124_prefix, n124_suffix, sep = "")

# n23
n23_prefix <- "bc"
n23_suffix <- seq(1:314278)
n23 = paste(n23_prefix, n23_suffix, sep = "")

# lambert
lambert <- c(n1234, n124)

# marioni
marioni <- c(n1234, n124, n23)

# jansen
n3_prefix <- "c"
n3_suffix <- seq(1:86818)
n3 = paste(n3_prefix, n3_suffix, sep = "")
jansen <- c(n1234, n23, n3)

# kunkle
n4_prefix <- "d"
n4_suffix <- seq(1:20391)
n4 = paste(n4_prefix, n4_suffix, sep = "")
kunkle <- c(n1234, n124, n4)

# myV <- plotVenn(list('Lambert (n = 74,046)' = lambert, 
#                      'Kunkle (n = 94,437)' = kunkle, 
#                      'Jansen (n = 455,258)' = jansen, 
#                      'Marioni (n = 388,324)' =  marioni), 
#                 nCycles = 6000)
myV <- plotVenn(list('Lambert' = lambert, 
                     'Kunkle' = kunkle, 
                     'Jansen' = jansen, 
                     'Marioni' =  marioni), 
                nCycles = 6000)
showSVG(myV, outFile = '~/Downloads/sampleoverlap.svg', 
        opacity = 0.1, borderWidth = 3, labelRegions = F, fontScale = 2.2, 
        setColors=c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A"))
fixsvg('~/Downloads/sampleoverlap.svg')

sprintf("inkscape -D -z -d 300 -w 744 --convert-dpi-method=scale-document --file=$PWD/%s --export-png=$PWD/%s", '~/Downloads/sampleoverlap.svg', '~/Downloads/sampleoverlap.png')

































