library(tidyverse)
library(ggstatsplot)
library(ggpubr)

setwd("/Users/Desktop/Dissertation/chapter_three/figures/in_vivo_splicing/")
cpm <- read.table(file = "cpm_df.txt", sep = "\t",
                  header = TRUE)




mt_one_v_two <- ggstatsplot::ggscatterstats(data = cpm, x = MT_Rep_one,
                                    y = MT_Rep_two,
                                    ) +
  labs(x = "MT_Rep_One",
       y = "MT_Rep_Two")
mt_two_v_three <- ggstatsplot::ggscatterstats(data = cpm, x = MT_Rep_two,
                                            y = MT_Rep_three,
) +
  labs(x = "MT_Rep_Two",
       y = "MT_Rep_Three")

pdf(file = "mt_cpm_correlation.pdf",
    width = 12, height = 6, paper = "letter")
ggarrange(mt_one_v_two, mt_two_v_three, ncol = 2, nrow = 1)
dev.off()


wt_one_v_two <- ggstatsplot::ggscatterstats(data = cpm, x = WT_Rep_one,
                            y = WT_Rep_two,
) +
  labs(x = "WT_Rep_one",
       y = "WT_Rep_two")
wt_two_v_three <- ggstatsplot::ggscatterstats(data = cpm, x = WT_Rep_two,
                                          y = WT_Rep_three,
) +
  labs(x = "WT_Rep_two",
       y = "WT_Rep_three")

pdf(file = "wt_cpm_correlation.pdf",
    width = 12, height = 6, paper = "letter")
ggarrange(wt_one_v_two, wt_two_v_three, ncol = 2, nrow = 1)
dev.off()
