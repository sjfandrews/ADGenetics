## Other tests 
pos <- position_jitter(width = 0.3, seed = 33)
ggplot(dat.p, aes(x = FRQ, y = neff, label = GENE, size = OR, fill = dir)) + 
  # geom_jitter(position = pos) +
  ggrepel::geom_label_repel(max.overlaps = 20, position = pos, size =2, min.segment.length = 10) +
  theme_bw()

ggsave("~/Dropbox/Research/PostDoc-MSSM/Neurogenomics/plots/neff_freq_labs.png")

pos <- position_jitter(width = 0.3, seed = 33)
ggplot(dat.p, aes(x = FRQ, y = dir, label = GENE)) + 
  # geom_jitter(position = pos) +
  ggrepel::geom_label_repel(max.overlaps = 20, position = pos, size =2, min.segment.length = 10) +
  theme_bw()

ggplot(dat.p, aes(x = FRQ, y = dir, label = GENE, size = dir, color = effect)) + 
  geom_point() +
  theme_bw() + scale_x_continuous(trans='log') + annotation_logticks(sides = "b")


ggplot(dat.p, aes(x = FRQ, y = dir, label = GENE, size = dir, color = effect)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(vars(effect), vars(maf), switch = "both", scales = "free") + 
  theme()


## Plot LD blocks

for(i in 1:nrow(ld)){
  message('Plotting ld for chr ', ld %>% slice(i) %>% pull(CHR), "...")
  
  if(!is.null(ld$ld[[i]])){
    df <- column_to_rownames(ld$ld[[i]], "RS_number") 
    rm_row = !rowSums(is.na(df)) == nrow(df)
    rm_col = !colSums(is.na(df)) == ncol(df)
    df <- df[rm_row, rm_col]
    
    p <- ggcorrplot::ggcorrplot(df, hc.order = TRUE, outline.col = "white")
    ggsave(glue("~/Dropbox/Research/PostDoc-MSSM/Neurogenomics/plots/ld_chr{chr}.png", chr = i))
    
  }
  
}
## Other Stuff 
adad.p <- ggplot(adad, aes(x = FRQ, y = OR, label = Gene)) + 
  geom_label_repel(min.segment.length = Inf, fill = "#D53E4F", size = geom.text.size) +
  # scale_x_continuous(limits = c(-0.1, 0.1), breaks=c(0), labels=c("Mutation")) +
  scale_y_continuous(limits = c(0, 1), breaks = 0.75, labels = "Causal") +
  theme_classic() + 
  labs(x = "Minor Allele Frequency", y = " ") + 
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.line.x=element_blank(), 
    axis.text.y = element_text(angle=90), 
    axis.ticks.y = element_blank(),
    text = element_text(size = theme.size)  
  )

adad.p


#### Colour by architecture
gwas.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  # geom_label(color = "black", show.legend = F) +
  # geom_text_repel(data = dat.p, aes(x = global_maf, y = OR, label = rare, point.size = size), 
  #                 segment.size  = 0.2, segment.color = "grey50", 
  #                 color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "risk" & GENE %nin% c("APOE","TREM2")), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(0.18, NA), 
                  segment.size  = 0.2, segment.color = "grey50", 
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "risk"), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(0.18, NA), 
                  segment.size  = 0.2, segment.color = "grey50", 
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "protective"), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(NA, -0.1), 
                  segment.size  = 0.2, segment.color = "grey50", 
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_point(data = dat.p %>% filter(architecture %in% c("Low frequency, Small effect", "Common, Small effect")),
             aes(x = global_maf, y = OR, color = architecture), size = 0.75) +
  geom_point(data = dat.p %>% filter(architecture %in% c("Low frequency, Moderate effect", "Common, Large effect")),
             aes(x = global_maf, y = OR, color = architecture, size = size)) +
  # geom_point(aes(size = size)) +
  scale_size(guide = 'none') +
  theme_classic() + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0005", "0.001", "0.01", "0.1", "0.25", "0.5")) + 
  scale_y_continuous(trans='log', limits = c(0.1, 4)) + 
  labs(x = "Global Population Minor Alelle Frequency", y = "Odds Ratio - Minor Allele") + 
  scale_color_manual(values = c("#377EB8","#4DAF4A",  "#FF7F00", "#984EA3")) +
  scale_fill_manual(values = c("#377EB8", "#4DAF4A", "#FF7F00", "#984EA3")) +
  guides(color=guide_legend(ncol=2)) + 
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.line.y=element_blank(), 
    text = element_text(size = theme.size)
  )

dat.p <- ad_loci %>% 
  group_by(locus) %>%
  arrange(year) %>%
  slice(which.max(year)) %>%
  ungroup() %>% 
  mutate(dir = ifelse(OR < 1, 1/OR, OR), 
         effect = case_when(
           OR > 3 ~ "large",
           between(OR, 1.3, 2.9) ~ "moderate", 
           OR < 1.3 ~ "small"
         ),
         effect = fct_relevel(effect, 'large', 'moderate', 'small'),
         maf = case_when(
           between(FRQ, 0.05, 1) ~ "common",
           between(FRQ, 0.001, 0.05) ~ "rare", 
           FRQ < 0.001 ~ "mutation"
         ), 
         maf = fct_relevel(maf, 'mutation', 'rare', 'common')
  ) %>%
  filter(!is.na(dir)) 


ggplot() + 
  # geom_jitter(data = filter(dat.p, dir == "risk") %>% rename(risk = OR), 
  #             aes(x = as_factor(year), y = FRQ, size = risk), 
  #             width = 0.25, colour = 'red') + 
  # geom_jitter(data = filter(dat.p, dir == "protective") %>% rename(protective = OR), 
  #             aes(x = as_factor(year), y = FRQ, size = protective), 
  #             width = 0.25, colour = 'blue') + 
  ggrepel::geom_label_repel(data = dat.p, aes(as_factor(year), y = FRQ, label = GENE), max.overlaps = 20) +
  theme_bw()

pos <- position_jitter(width = 0.3, seed = 33)
ggplot(data = dat.p, aes(as_factor(year), y = FRQ, label = GENE, size = OR, fill = dir)) + 
  # geom_jitter(position = pos) +
  ggrepel::geom_label_repel(max.overlaps = 20, position = pos, size =2, min.segment.length = 10) +
  theme_bw()

ggsave("~/Dropbox/Research/PostDoc-MSSM/Neurogenomics/plots/freq_year_lab.png")



### Ternary color scheme
# Ternary
# gg = ggplotGrob(tric_educ$key + theme_showarrows() + theme_hidelabels() + theme_notitles() + 
#                   theme(text = element_text(size = theme.size))) 
# 
# plotGG(
#   plot = gg,
#   x = 9, y = -0.1,
#   width = 1.7, height = 1.7, just = c("right", "top")
# )

ADAD <- tribble(
  ~GENE, ~OR, ~global_maf, ~annotation_impact,
  "ADAD", 10, 0, "PATHOGENIC",
)

std_100 <- function(x){(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))}

(0.06315 - 0) / (0.49 - 0)

test <- dat.p %>% 
  select(GENE, OR, global_maf, annotation_impact, architecture) %>%
  bind_rows(ADAD) %>%
  mutate(
    OR2 = ifelse(OR < 1, 1/OR, OR),
    effect = std_100(OR2),
    freq = std_100(global_maf),
    effect = effect * 100, 
    freq = freq * 100,
    impact = case_when(
      architecture  == "Common, Small effect" ~ 0,
      annotation_impact == "PATHOGENIC" ~ 100,
      annotation_impact == "HIGH" ~ 75,
      annotation_impact == "MODERATE" ~ 50,
      annotation_impact == "LOW" ~ 25,
      annotation_impact == "MODIFIER" ~ 0,
    ), 
    # effect = (abs(impact) - min(abs(impact), na.rm = T)) / (max(abs(impact), na.rm = T) - min(abs(impact), na.rm = T)),
  )

tri_center = rep(1/3, 3)
# tri_center = c(0.49, 0.02, 0.49)
# tri_center = c(0.01, 0.01, 0.98)

tric_educ <- Tricolore(test,
                       p1 = 'effect', p2 = 'freq', p3 = 'impact', 
                       breaks = Inf, crop = FALSE, legend = TRUE, spread = )

tric_educ$key + theme_showarrows() + theme_hidelabels() + theme_notitles() + geom_crosshair_tern()

df = data.frame(x=100,y=0,z=100)
base = ggtern(df,aes(x,y,z)) + geom_point()
base + geom_crosshair_tern()

test$tri <- tric_educ$rgb
test$names <- as.character(1:nrow(test))
test_colours <- setNames(test$tri, test$names)

dat.p$tri <- tric_educ$rgb[-length(tric_educ$rgb)] 
dat.p$names <- as.character(1:nrow(dat.p))

colours <- setNames(dat.p$tri, dat.p$names)

gwas.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  # geom_label(color = "black", show.legend = F) +
  # geom_text_repel(data = dat.p, aes(x = global_maf, y = OR, label = rare, point.size = size), 
  #                 segment.size  = 0.2, segment.color = "grey50", 
  #                 color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "risk"), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(0.18, NA), 
                  segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "protective"), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(NA, -0.1), 
                  segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_point(data = dat.p %>% filter(architecture %in% c("Low frequency, Small effect", "Common, Small effect")),
             aes(x = global_maf, y = OR, color = names), size = 0.75) +
  geom_point(data = dat.p %>% filter(architecture %in% c("Low frequency, Moderate effect", "Common, Large effect")),
             aes(x = global_maf, y = OR, color = names, size = size)) +
  # geom_point(aes(size = size)) +
  scale_size(guide = 'none') +
  theme_classic() + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0005", "0.001", "0.01", "0.1", "0.25", "0.5")) + 
  scale_y_continuous(trans='log', limits = c(0.1, 8)) + 
  labs(x = "Global Population Minor Alelle Frequency", y = "Odds Ratio - Minor Allele") + 
  scale_colour_manual(values = colours) + 
  guides(color=guide_legend(ncol=2)) + 
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.line.y=element_blank(), 
    text = element_text(size = theme.size), 
    legend.position = 'none'
  )

gwas.p 


#### color by vep
# Normal
# plotGG(
#   plot = out.legend,
#   x = 1.5, y = 0,
#   width = 7.5, height = 0.75, just = c("left", "top")
# )


gwas.p <- ggplot() + 
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  # geom_label(color = "black", show.legend = F) +
  # geom_text_repel(data = dat.p, aes(x = global_maf, y = OR, label = rare, point.size = size), 
  #                 segment.size  = 0.2, segment.color = "grey50", 
  #                 color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "risk"), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(0.18, NA), 
                  segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_text_repel(data = filter(dat.p, dir == "protective"), 
                  aes(x = global_maf, y = OR, label = label, point.size = size), 
                  max.overlaps = 20, ylim = c(NA, -0.1), 
                  segment.size  = 0.2, segment.color = "grey50", min.segment.length = 0,
                  color = "black", size = geom.text.size, show.legend = F) +
  geom_point(data = dat.p %>% filter(architecture %in% c("Low frequency, Small effect", "Common, Small effect")),
             aes(x = global_maf, y = OR, color = global_maf), size = 0.75) +
  geom_point(data = dat.p %>% filter(architecture %in% c("Low frequency, Moderate effect", "Common, Large effect")),
             aes(x = global_maf, y = OR, color = global_maf, size = size)) +
  # geom_point(aes(size = size)) +
  scale_size(guide = 'none') +
  theme_classic() + 
  scale_x_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5),
                     labels = c("0.0005", "0.001", "0.01", "0.1", "0.25", "0.5")) + 
  scale_y_continuous(trans='log', limits = c(0.1, 8)) + 
  labs(x = "Global Population Minor Alelle Frequency", y = "Odds Ratio - Minor Allele") + 
  # scale_colour_gradient(low = "#D53E4F", high = "#3288BD", trans = "log") +
  # scale_colour_gradient2(low = "#D53E4F", mid = "#762a83", high = "#3288BD", midpoint = 0.05, trans = "log") +
  # colorspace::scale_color_continuous_sequential(palette = 'plasma') +
  scale_color_continuous_divergingx(palette = "Zissou1", mid = -5, rev = T, trans = 'log') +
  # scale_color_manual(values = c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")) +
  # scale_fill_manual(values = c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")) +
  guides(color=guide_legend(ncol=2)) + 
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.line.y=element_blank(), 
    text = element_text(size = theme.size)
  )

gwas.p 
