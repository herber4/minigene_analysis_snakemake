library(tidyverse)
library(ggpubr)
library(wordcloud2)
library(ggrain)
library(ggbeeswarm)
setwd("/Users/herber4/Desktop/Dissertation/chapter_three/figures/in_vivo_splicing/")

dpsi <- read.table(file = "/Users/herber4/Desktop/Dissertation/Dissertation_Aggregation/github_data/Figure_4_data.txt",
                   sep = "\t", header = TRUE)

tmp <- read.table(file = "K700E_vs_WT_dPSI_fully_annotated.txt",
                   sep = "\t", header = TRUE)
tmp <- tmp[,c(1,4,5)]
dpsi$dpsi_norm <- dpsi$delta_PSI-.105
dpsi <- dpsi %>%
  mutate(color = dpsi_norm >= 0.05)
dpsi <- merge(dpsi, tmp, by = "sample")
pdf(file = "../figure_5/K700E_vs_WT_BPS_SimScore_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
dpsi %>%
  filter(str_detect(annot, "Branch_Point")) %>%
  ggplot(aes(x = delta_PSI, y = similarity_score, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(x = "dC3SS (K700E - WT)",
       y = "Similarity Score (RNAsmc)") +
  theme_bw() +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()

model <- lm(similarity_score ~ delta_PSI, data = dpsi %>% filter(str_detect(annot, "Branch_Point")))

summary(model)

pdf(file = "../figure_5/K700E_vs_WT_RBP_Block_Three_SimScore_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
dpsi %>%
  filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three")) %>%
  ggplot(aes(x = delta_PSI, y = similarity_score, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (K700E - WT)",
       y = "Similarity Score (RNAsmc)") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()

model <- lm(similarity_score ~ delta_PSI, data = dpsi %>% filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three")))

summary(model)

pdf(file = "../figure_5/K700E_vs_WT_RBP_Block_One_SimScore_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
dpsi %>%
  filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One")) %>%
  ggplot(aes(x = delta_PSI, y = similarity_score, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (K700E - WT)",
       y = "Similarity Score (RNAsmc)") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()
 model <- lm(similarity_score ~ delta_PSI, data = dpsi %>% filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One")))

summary(model)


wtvsmt <- read.table(file = "../../../Dissertation_Aggregation/github_data/MTRNA_vs_WTRNA_Statistics_table.txt",
                     sep = "\t", header = TRUE)

wtvsmt <- merge(wtvsmt, tmp, by = "sample")
tmp2 <- dpsi[,c(1,3)]
wtvsmt <- merge(wtvsmt, tmp2, by = "sample")
wtvsmt <- wtvsmt %>%
  mutate(signif = p_adj < .05)


pdf(file = "../figure_5/WT_vs_MT_RBP_Block_One_SimScore_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
wtvsmt %>%
  filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One")) %>%
  ggplot(aes(x = delta_C3SS_ratio, y = similarity_score, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (MT - WT)",
       y = "Similarity Score (RNAsmc)") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()

model <- lm(similarity_score ~ delta_C3SS_ratio, data = wtvsmt %>% filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One")))

summary(model)


pdf(file = "../figure_5/WT_vs_MT_RBP_Block_Three_SimScore_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
wtvsmt %>%
  filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three")) %>%
  ggplot(aes(x = delta_C3SS_ratio, y = similarity_score, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (MT - WT)",
       y = "Similarity Score (RNAsmc)") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()

model <- lm(similarity_score ~ delta_C3SS_ratio, data = wtvsmt %>% filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three")))

summary(model)



pdf(file = "../figure_5/WT_vs_MT_Branch_Point_SimScore_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
wtvsmt %>%
  filter(str_detect(annot, "Branch_Point")) %>%
  ggplot(aes(x = delta_C3SS_ratio, y = similarity_score, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (MT - WT)",
       y = "Similarity Score (RNAsmc)") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()

model <- lm(similarity_score ~ delta_C3SS_ratio, data = wtvsmt %>% filter(str_detect(annot, "Branch_Point")))

summary(model)
