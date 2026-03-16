library(tidyverse)
library(ggpubr)
library(wordcloud2)
library(ggrain)
library(ggbeeswarm)
setwd("/Users/herber4/Desktop/Dissertation/chapter_three/figures/in_vivo_splicing/")

dpsi <- read.table(file = "K700E_vs_WT_dPSI_fully_annotated.txt",
                  sep = "\t", header = TRUE)

dpsi$dpsi_norm <- dpsi$delta_PSI-.105
dpsi <- dpsi %>%
  mutate(color = dpsi_norm >= 0.05)

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


### now for SHAPE correlations

pdf(file = "../figure_5/K700E_vs_WT_BPS_SHAPE_COR_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
dpsi %>%
  filter(str_detect(annot, "Branch_Point")) %>%
  ggplot(aes(x = delta_PSI, y = spearman_cor, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(x = "dC3SS (MT oligo - WT oligo)",
       y = "Spearman_Corr") +
  theme_bw() +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()

pdf(file = "../figure_5/K700E_vs_WT_RBP_Block_Three_SHAPE_COR_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
dpsi %>%
  filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three")) %>%
  ggplot(aes(x = delta_PSI, y = spearman_cor, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (MT oligo - WT oligo)",
       y = "Spearman_Corr") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()


pdf(file = "../figure_5/K700E_vs_WT_RBP_Block_One_SHAPE_SCORE_vs_dC3SS.pdf",
    height = 12, width = 10, paper = "letter")
dpsi %>%
  filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One")) %>%
  ggplot(aes(x = delta_PSI, y = spearman_cor, col = signif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(x = "dC3SS (MT oligo - WT oligo)",
       y = "Spearman_Corr") +
  scale_color_manual(values = c("TRUE"="blue", "FALSE"="red"))
dev.off()



### fetch sequences

bps <- dpsi %>%
  filter(str_detect(annot, "Branch_Point"))

bps$source <- "Branch_Point"

rbp_three <- dpsi %>%
  filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three"))
rbp_three$source <- "RBP_Three"

rbp_one <- dpsi %>%
  filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One"))
rbp_one$source <- "RBP_One"

fetch <- rbind(rbp_one, bps, rbp_three)

seqs <- read.table(file = "../figure_5/all_seqs.txt",
                   sep = "\t", header = TRUE)

seqs <- merge(fetch, seqs, by = "sample")
write.table(seqs, file = "../figure_5/K700E_vs_WT_RBP_One_Three_and_BPS_Sequences.txt",
            sep = "\t", row.names = F, quote = F)

