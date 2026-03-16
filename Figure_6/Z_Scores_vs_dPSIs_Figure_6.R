library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(factoextra)
library(mclust)
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

wt_scan <- merge(wt_scan, annots, by = c("sample", "Nucleotide"))
wt_scan[,c(19:39)] <- NULL

tmp <- rbind(wt_scan, tmp)


red <- tmp %>%
  group_by(source) %>%
  summarise(mean_z = mean(avgZ),
            mean_ED = mean(avgED),
            z_sd = sd(avgZ),
            ed_sd = sd(avgED),
            k700e_dpsi = mean(K700E_vs_WT_dPSI),
            wt_vs_mt_dpsi = mean(MT_vs_WT_dPSI),
            similarity_score = mean(similarity_score),
            correlation = mean(spearman_cor, na.rm = TRUE))

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


red2 %>%
  filter(!source == "NA") %>%
  ggplot(aes(x = k700e_dpsi, y = mean_z, col = source)) +
  geom_point()

red2 %>%
  filter(!source == "NA") %>%
  ggplot(aes(x = k700e_dpsi, y = mean_z, col = source)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ source)


library(ggplot2)

red2 %>%
  filter(!source == "NA",
         !source == "WT") %>%
  ggplot(aes(x = k700e_dpsi, y = mean_z, color = source)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "K700E dPSI",
    y = "Mean Z Score",
    color = "Source"
  )

red2 %>%
  filter(!source == "NA",
         !source == "WT") %>%
  ggplot(aes(x = wt_vs_mt_dpsi, y = mean_z, color = source)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "MT_vs_WT_dPSI",
    y = "Mean Z Score",
    color = "Source"
  )



ggplot(red2, aes(x = mean_z, y = k700e_dpsi)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ source) +
  theme_classic() +
  labs(
    x = "Mean Z",
    y = "K700E dPSI"
  )


model <- lm(k700e_dpsi ~ mean_z * source, data = red2)
summary(model)

ggplot(red, aes(x = wt_vs_mt_dpsi, y = mean_z, col = source)) +
  geom_point(size = 2) +
  theme_bw()


dpsi <- read.table(file = "../in_vivo_splicing/K700E_vs_WT_dPSI_fully_annotated.txt",
                   sep = "\t", header = TRUE)

ref <- dpsi %>%
  filter(sample == "Two_205")
dpsi <- dpsi %>%
  filter(str_detect(annot, "Branch_Point|Stem_4_2|RBP_Block_Three|RBP_three|Stem_3_1|RBP_Block_One|RBP_One"))


tmp <- tmp[,c(1:8, 24)]
wt_scan$source <- "WT"

tmp <- rbind(tmp, wt_scan)
