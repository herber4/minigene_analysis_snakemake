library(tidyverse)
library(ggpubr)
setwd("/Users/herber4/Desktop/Dissertation/chapter_three/figures/in_vivo_splicing/")

dpsi <- read.table(file = "../../../Dissertation_Aggregation/github_data/K700ERNA_vs_WTRNA_Stats_Table.txt",
                   sep = "\t", header = TRUE)
dpsi <- dpsi %>%
  mutate(
    signif = padj < 0.05
  )
dpsi$delta_PSI <- dpsi$MT-dpsi$WT
pdf(file = "dPSI_K700E_vs_WT_per_Oligo_Volcanoe.pdf",
    height = 12, width = 12, paper = "letter")
ggplot(dpsi, aes(x = delta_PSI, y = -log10(padj), color = signif)) +
  geom_point() +
  geom_text(
    data = subset(dpsi, sample == "Two_189"),
    aes(label = sample),
    vjust = -0.5
  ) +
  scale_color_manual(
    values = c(`TRUE` = "blue", `FALSE` = "black")
  ) +
  theme_bw() +
  labs(
    x = expression(Delta*C3SS~"(MT - WT)"),
    y = expression(-log[10]("FDR")),
    color = "FDR < 0.05"
  ) +
  geom_vline(xintercept = .09, linetype = "dashed", col = "red") +
  theme(legend.position = "")
dev.off()

hist(dpsi$delta_PSI)


#### plot dpsi versus RNAsmc structure similarity

sims <- read.table(file = "structure_similarity_scores.txt",
                   sep = "\t", header = TRUE)

sims$samps <- gsub(".ct", "", sims$samps)
sims$samples <- gsub(".ct", "", sims$samples)

sims %>%
  filter(samps == "One_205") %>%
  ggplot(aes(x = reorder(samples, -similarity_score), y = similarity_score)) +
  geom_col()

ref <- sims %>%
  filter(samps == "One_205")
colnames(ref)[2] <- "sample"
ref$samps <- NULL

dpsi <- merge(ref, dpsi, by = "sample")

pdf(file = "K700E_vs_WT_SimScore_vs_DPSI.pdf",
    height = 12, width = 12, paper = "letter")
dpsi %>%
  filter(signif == TRUE) %>%
  ggplot(aes(x = delta_PSI, y = similarity_score)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw()
dev.off()
model <- lm(similarity_score ~ delta_PSI, data = dpsi %>% filter(signif == "TRUE"))
summary(model)
### now do spearman correlation of SHAPE vs One_205

shape <- read.table(file = "/Users/herber4/Desktop/Dissertation/chapter_three/data/NIA_Pilot_vs_Pooled_Shape_Data_For_Figures.txt",
                    sep = "\t", header = TRUE)
hist(dpsi$delta_PSI)
shape <- shape %>%
  filter(rep == "One")
shape$rep <- NULL
colnames(shape)[30] <- "sample"

s_ref <- shape %>%
  filter(sample == "One_205")
shape <- shape %>%
  filter(!sample == "One_205")
library(dplyr)

library(dplyr)

shape_cor <- shape %>%
  left_join(
    s_ref %>%
      select(Nucleotide, ref_profile = Norm_profile),
    by = "Nucleotide"
  ) %>%
  group_by(sample) %>%
  summarise(
    spearman_cor = {
      df <- cur_data()
      if (sum(complete.cases(df$Norm_profile, df$ref_profile)) >= 2) {
        cor(
          df$Norm_profile,
          df$ref_profile,
          method = "spearman",
          use = "complete.obs"
        )
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

dpsi <- merge(shape_cor, dpsi, by = "sample")

pdf(file = "K700E_vs_WT_SHAPE_Corr_vs_dPSI.pdf",
    width = 12, height = 12, paper = "letter")
dpsi %>%
  filter(signif == TRUE) %>%
  ggplot(aes(x = delta_PSI, y = spearman_cor)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw()
dev.off()

model <- lm(spearman_cor ~ delta_PSI, data = dpsi %>% filter(signif == "TRUE"))
summary(model)
