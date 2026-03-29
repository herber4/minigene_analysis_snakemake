library(tidyverse)

setwd(
  "/Users/herber4/Desktop/Dissertation/chapter_three/figures/in_vivo_splicing/"
)

dpsi <- read.table(file = "dPSI_WToligo_vs_MToligo.txt",
                   sep = "\t", header = TRUE)


dpsi <- dpsi %>%
  mutate(
    signif = p_adj < 0.05
  )


annots <- read.table(file = "dPSI_Annotated.txt",
                   sep = "\t", header = TRUE)

annots <- annots[,c(1:3)]

dpsi <- merge(annots, dpsi, by = "sample")


library(ggplot2)

pdf(file = "deltaPSI_WToligo_vs_MToligo_Volcano.pdf",
    width = 12, height = 8, paper = "letter")
ggplot(dpsi, aes(x = delta_C3SS_ratio, y = -log10(p_adj), color = signif)) +
  geom_point() +
  geom_text(aes(label = ifelse(dpsi$name == "Rep_Two_MAP_Can_CAG_to_ACG", as.character(dpsi$annot), ""),
                hjust = 0, vjust = 0)) +
  geom_text(aes(label = ifelse(dpsi$name == "Rep_Two_MAP_Weaken_Canonical_PyTract", as.character(dpsi$annot), ""),
                hjust = 0, vjust = 0)) +
  scale_color_manual(
    values = c(`TRUE` = "blue", `FALSE` = "black")
  ) +
  theme_bw() +
  labs(
    x = expression(Delta*C3SS~"(MT Oligo - WT Oligo)"),
    y = expression(-log[10]("FDR")),
    color = "p value < 0.05"
  ) 
dev.off()

### this is structure similarity scores vs delta PSI

sims <- read.table(file = "structure_similarity_scores.txt",
                   sep = "\t", header = TRUE)

sims$samps <- gsub(".ct", "", sims$samps)
sims$samples <- gsub(".ct", "", sims$samples)

sims <- sims %>%
  filter(samps == "One_205")
colnames(sims)[2] <- "sample"
sims$samps <- NULL 

dpsi <- merge(dpsi, sims, by = "sample")

pdf(file = "WToligo_vs_MToligo_RNAsmc_over_dPSI.pdf",
    width = 12, height = 8, paper = "letter")
dpsi %>%
  filter(signif == "TRUE") %>%
  ggplot(aes(x = delta_C3SS_ratio, y = similarity_score)) +
  geom_point() +
  theme_bw() +
  geom_smooth(data = dpsi %>% filter(signif == "TRUE"), method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(
    x = expression(Delta*C3SS~"(MT Oligo - WT Oligo)"),
    y = "Structure Similarity Score (RNAsmc)"
  ) 
dev.off()

model <- lm(similarity_score ~ delta_C3SS_ratio, data = dpsi %>% filter(signif == "TRUE"))

summary(model)

### this is SHAPE correlation to deltaPSI

shape <- read.table(file = "/Users/herber4/Desktop/Dissertation/chapter_three/data/NIA_Pilot_vs_Pooled_Shape_Data_For_Figures.txt",
                    sep = "\t", header = TRUE)
shape <- shape %>%
  filter(rep == "One")
shape <- shape[,c(1, 28, 30)]
colnames(shape)[3] <- "sample"

library(dplyr)

ref <- shape %>%
  filter(sample == "One_205") %>%
  select(Nucleotide, ref_profile = Norm_profile)

shape_cor <- shape %>%
  filter(sample != "One_205") %>%
  left_join(ref, by = "Nucleotide") %>%
  group_by(sample) %>%
  summarise(
    n_used = sum(complete.cases(Norm_profile, ref_profile)),
    spearman_cor = ifelse(
      n_used >= 3,  # minimum reasonable number
      cor(
        Norm_profile,
        ref_profile,
        method = "spearman",
        use = "complete.obs"
      ),
      NA_real_
    ),
    .groups = "drop"
  )

shape_cor <- bind_rows(
  tibble(sample = "One_205", spearman_cor = 1),
  shape_cor
)


dpsi <- merge(dpsi, shape_cor, by = "sample")

pdf(file = "WToligo_vs_MToligo_SHAPEcorr_over_dPSI.pdf",
    width = 12, height = 8, paper = "letter")
dpsi %>%
  filter(signif == "TRUE") %>%
  ggplot(aes(x = delta_C3SS_ratio, y = spearman_cor)) +
  geom_point() +
  theme_bw() +
  geom_smooth(data = dpsi %>% filter(signif == "TRUE"),method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(x = expression(Delta*C3SS~"(MT Oligo - WT Oligo)"),
       y = "Spearman Correlation (SHAPE)")
dev.off()

model <- lm(spearman_cor ~ delta_C3SS_ratio, data = dpsi %>% filter(signif == "TRUE"))

summary(model)


pdf(file = "WToligo_vs_MToligo_SHAPE_Corr_vs_RNAsmc.pdf",
    width = 12, height = 8, paper = "letter")
dpsi %>%
  filter(signif == "TRUE") %>%
  ggplot(aes(x = similarity_score, y = spearman_cor)) +
  geom_point() +
  geom_smooth(data = dpsi %>% filter(signif == "TRUE"),method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(x = "Similarity Score (RNAsmc)",
       y = "Spearman Correlation (SHAPE)") +
  theme_bw()
dev.off()

model <- lm(spearman_cor ~ similarity_score, data = dpsi %>% filter(signif == "TRUE"))

summary(model)
