library(dplyr)
library(RColorBrewer)
library(circlize)
library(ggpubr)
library(stringr)
library(ggrain)
setwd("/Users/herber4/Desktop/Dissertation/chapter_three/figures/scanfold/")

annots <- read.table(file = "../../claude/data_4_claude.txt",
                     sep = "\t", header = TRUE)
scan <- read.table(file = "scanfold_per_base_Z_scores.txt",
                   sep = "\t", header = TRUE)

one <- read.table(file = "One_205_stats.txt",
                  sep = "\t", header = TRUE)
one$source <- NULL

scan$gene <- sub("\\..*", "", scan$gene)
colnames(scan)[7] <- "sample"

annots <- rbind(annots, one)

mas <- merge(scan, annots,  by = c("sample", "Nucleotide"))

mas$k700e_signif <- ifelse(mas$K700E_vs_WT_p_value < .05, "True", "False")
mas$WT_signif <- ifelse(mas$MT_vs_WT_p_value < .05, "True", "False")


wt_scan <- scan %>%
  filter(str_detect(sample, "205"))

wt_scan <- wt_scan %>%
  filter(Nucleotide >= 56, Nucleotide <= 225)

mas[,c(18:38)] <- NULL

bps <- mas %>%
  filter(str_detect(annot, "Branch_Point"))

bps$source <- "Branch_Point"

rbp_three <- mas %>%
  filter(str_detect(annot, "Stem_4_2|RBP_Block_Three|RBP_three"))
rbp_three$source <- "RBP_Three"

rbp_one <- mas %>%
  filter(str_detect(annot, "Stem_3_1|RBP_Block_One|RBP_One"))
rbp_one$source <- "RBP_One"

wt <- mas %>%
  filter(str_detect(annot, "WT"))

wt$source <- "WT"

tmp <- rbind(rbp_one, bps, rbp_three, wt)

tmp2 <- anti_join(mas, tmp, by = "sample")

tmp2$source <- "NA"

tmp <- rbind(tmp, tmp2)

rm(tmp2, wt, rbp_one, rbp_three, bps)

#tmp <- tmp[,c(1:8, 24)]

wt_scan$source <- "WT"
wt_scan$k700e_signif <- "True"
wt_scan$WT_signif <- "True"
wt_scan <- merge(wt_scan, annots, by = c("sample", "Nucleotide"))


wt_scan[,c(21:41)] <- NULL

tmp <- rbind(wt_scan, tmp)

tmp$WT_signif[is.na(tmp$WT_signif)] <- "True"


k_stats <- tmp %>%
  group_by(source, k700e_signif, Nucleotide) %>%
  summarise(mean_z = mean(avgZ),
            mean_ED = mean(avgED),
            z_sd = sd(avgZ),
            ed_sd = sd(avgED),
            k700e_dpsi = mean(K700E_vs_WT_dPSI),
            wt_vs_mt_dpsi = mean(MT_vs_WT_dPSI),
            similarity_score = mean(similarity_score),
            correlation = mean(spearman_cor, na.rm = TRUE))

wt_stats <- tmp %>%
  group_by(source, WT_signif, Nucleotide) %>%
  summarise(mean_z = mean(avgZ),
            mean_ED = mean(avgED),
            z_sd = sd(avgZ),
            ed_sd = sd(avgED),
            k700e_dpsi = mean(K700E_vs_WT_dPSI),
            wt_vs_mt_dpsi = mean(MT_vs_WT_dPSI),
            similarity_score = mean(similarity_score),
            correlation = mean(spearman_cor, na.rm = TRUE))
dpsi_wt <- tmp %>%
  group_by(source, WT_signif) %>%
  summarise(k700e_dpsi = mean(K700E_vs_WT_dPSI),
            wt_vs_mt_dpsi = mean(MT_vs_WT_dPSI),
            similarity_score = mean(similarity_score),
            correlation = mean(spearman_cor, na.rm = TRUE))


dpsi_k <- tmp %>%
  group_by(source, k700e_signif) %>%
  summarise(k700e_dpsi = mean(K700E_vs_WT_dPSI),
            wt_vs_mt_dpsi = mean(MT_vs_WT_dPSI),
            similarity_score = mean(similarity_score),
            correlation = mean(spearman_cor, na.rm = TRUE))


k_stats %>%
  filter(source == "RBP_Three") %>%
  ggplot(aes(x = factor(Nucleotide), y = mean_z, col = k700e_signif, group = k700e_signif)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = factor(189), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = factor(209), linetype = "dashed", color = "black") +
  theme_bw()


#MT oligo vs WT oligo

tmp$group_fill_w <- interaction(tmp$source, tmp$WT_signif)


pdf(file = "MT_vs_WT_oligo_Signif_vs_Non_Z_score_boxes.pdf")


box <- ggplot(tmp, aes(x = source, y = avgZ, fill = group_fill_w)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group",
       y = "Z Score",
       fill = "Significant") +
  scale_fill_manual(values = c(
    "RBP_One.False"      = "#9ecae1",
    "RBP_One.True"       = "#08519c",
    "RBP_Three.False"    = "#a1d99b",
    "RBP_Three.True"     = "#006d2c",
    "Branch_Point.False" = "#fcae91",
    "Branch_Point.True"  = "#cb181d",
    "WT.True"            = "#54278f",
    "NA.False"           = "#d9d9d9",
    "NA.True"            = "#636363"
  )) +
  theme(legend.position = "")

pdf(file = "new_box_legend.pdf",
    width = 12, height = 4, paper = "letter")
ggplot(tmp, aes(x = source, y = avgZ, fill = group_fill_w)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group",
       y = "Z Score",
       fill = "Significant") +
  scale_fill_manual(values = c(
    "RBP_One.False"      = "#9ecae1",
    "RBP_One.True"       = "#08519c",
    "RBP_Three.False"    = "#a1d99b",
    "RBP_Three.True"     = "#006d2c",
    "Branch_Point.False" = "#fcae91",
    "Branch_Point.True"  = "#cb181d",
    "WT.True"            = "#54278f",
    "NA.False"           = "#d9d9d9",
    "NA.True"            = "#636363"
  ))
dev.off()

red2 <- tmp %>%
  group_by(source, sample) %>%
  summarise(mean_z = mean(avgZ),
            mean_ED = mean(avgED),
            z_sd = sd(avgZ),
            ed_sd = sd(avgED),
            k700e_dpsi = mean(K700E_vs_WT_dPSI),
            wt_vs_mt_dpsi = mean(MT_vs_WT_dPSI),
            similarity_score = mean(similarity_score),
            correlation = mean(spearman_cor, na.rm = TRUE))


dot <- red2 %>%
  filter(!source == "NA",
         !source == "WT") %>%
  ggplot(aes(x = k700e_dpsi, y = mean_z, color = source)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "K700E dPSI",
    y = "Mean Z Score",
    color = "Source"
  ) +
  theme(legend.position = "")

pdf(file = "MT_vs_WT_Signif_vs_Non_Box_N_Dot_plot.pdf",
    width = 12, height = 6, paper = "letter")
ggarrange(box, dot, nrow = 1, ncol = 2)
dev.off()

pdf(file = "new_line_legend.pdf",
    width = 12, height = 4, paper = "letter")
red2 %>%
  filter(!source == "NA",
         !source == "WT") %>%
  ggplot(aes(x = k700e_dpsi, y = mean_z, color = source)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, ) +
  stat_cor() +
  theme_classic() +
  labs(
    x = "K700E dPSI",
    y = "Mean Z Score",
    color = "Source"
  )
dev.off()


pairwise.wilcox.test(tmp$avgZ,
                     tmp$group_fill_w,
                     p.adjust.method = "BH")
