library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(factoextra)
library(mclust)
library(RColorBrewer)
library(circlize)
library(ggpubr)
setwd("/Users/herber4/Desktop/Dissertation/chapter_three/figures/scanfold/")

annots <- read.table(file = "../../claude/data_4_claude.txt",
                 sep = "\t", header = TRUE)

scan <- read.table(file = "scanfold_per_base_Z_scores.txt",
                   sep = "\t", header = TRUE)



scan$gene <- sub("\\..*", "", scan$gene)
colnames(scan)[7] <- "sample"

mas <- merge(scan, annots,  by = c("sample", "Nucleotide"))

wt_scan <- scan %>%
  filter(sample == "One_205")

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


tmp <- tmp[,c(1:8, 24)]
wt_scan$source <- "WT"

tmp <- rbind(tmp, wt_scan)

stats <- tmp %>%
  group_by(source, Nucleotide) %>%
  summarise(mean_z = mean(avgZ),
            mean_ED = mean(avgED))


z_line <- ggplot(stats, aes(x = factor(Nucleotide), y = mean_z, col = source, group = source)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = factor(189), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = factor(209), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(x = "Position", y = "Mean Z Score") +
  theme(axis.text.x = element_blank(),
        legend.position = "")

ed_line <- ggplot(stats, aes(x = factor(Nucleotide), y = mean_ED, col = source, group = source)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = factor(189), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = factor(209), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(y = "Mean Ensemble Diversity",
       x = "") +
  theme(axis.text.x = element_blank(),
        legend.position = "")

ed_line

pdf(file = "liners.pdf",
    width = 16, height = 6, paper = "letter")
ggarrange(ed_line, z_line, ncol = 1, nrow = 2)
dev.off()

z_box <- ggplot(tmp, aes(x = source, y = avgZ, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group",
       y = "Z Score") +
  theme(legend.position = "")


ed_box <- ggplot(tmp, aes(x = source, y = avgED, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group",
       y = "Ensemble Diversity") +
  theme(legend.position = "")


ed_box

pdf(file = "Boxes.pdf",
    width = 12, height = 6, paper = "letter")
ggarrange(ed_box, z_box, ncol = 2, nrow = 1)
dev.off()

pdf(file = "box_legend.pdf",
    width = 12, height = 10, paper = "letter")
ggplot(tmp, aes(x = source, y = avgZ, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group",
       y = "Z Score")
dev.off()

pdf(file = "line_legend.pdf",
    width = 12, height = 6, paper = "letter")
ggplot(stats, aes(x = factor(Nucleotide), y = mean_ED, col = source, group = source)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = factor(189), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = factor(209), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(y = "Mean Ensemble Diversity",
       x = "")
dev.off()


pairwise.t.test(tmp$avgED, tmp$source,
                p.adjust.method = "BH")

pairwise.t.test(tmp$avgZ, tmp$source,
                p.adjust.method = "BH")

pairwise.t.test(stats$mean_z, stats$source,
                p.adjust.method = "BH")
