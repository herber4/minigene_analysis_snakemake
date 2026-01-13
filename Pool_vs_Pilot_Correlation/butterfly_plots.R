library(tidyverse)
library(ggpubr)

shape <- read.table(file = "/Desktop/Dissertation/chapter_three/data/NIA_Pilot_vs_Pooled_Shape_Data_For_Figures.txt",
                    sep = "\t", header = TRUE)
shape %>%
  filter(source %in% c("base_base", "One_205", "Two_205", "bp_bp", "disrupted_disrupted")) %>%
  ggplot(aes(x = as.factor(Nucleotide), y = Norm_profile, group = source, col = source)) +
  geom_line()


shape %>%
  filter(source %in% c("base_base", "One_205", "Two_205")) %>%
  ggplot(aes(x = as.factor(Nucleotide), y = Norm_profile, group = source, col = source)) +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept = as.factor(188), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = as.factor(208), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("base_base"="#e31a1c", "One_205"="#b3cde3",
                                "Two_205"="#c2e699")) +
  theme(axis.text.x = element_blank()) +
  labs(x = "Nucleotide")


pdf(file = "/Desktop/Dissertation/chapter_three/figures/NIA_Pilot_vs_Pool/WT_Pilot_vs_WT_Pooled.pdf",
    width = 12, height = 4, paper = "letter")
wt <- shape %>%
  filter(source == "base_base") %>%
  ggplot(aes(x = as.factor(Nucleotide), y = Norm_profile)) +
  geom_col(fill = "#f1b6da") +
  geom_col(
    data = shape %>% filter(source == "One_205"),
    aes(x = as.factor(Nucleotide), y = -Norm_profile),
    fill = "#b8e186"
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_vline(xintercept = as.factor(188), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = as.factor(208), color = "black", linetype = "dashed") +
  ggtitle("WT MAP3K7 Pilot vs Matched Pooled Oligo") +
  labs(x = "Nucleotide, Position",
       y = "Normalized SHAPE Reactivity")


#### pearson, vs Two_181 == .73
bps <- shape %>%
  filter(source == "bp_bp") %>%
  ggplot(aes(x = as.factor(Nucleotide), y = Norm_profile)) +
  geom_col(fill = "#f1b6da") +
  geom_col(
    data = shape %>% filter(source == "Two_181"),
    aes(x = as.factor(Nucleotide), y = -Norm_profile),
    fill = "#b8e186"
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_vline(xintercept = as.factor(188), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = as.factor(208), color = "black", linetype = "dashed") +
  ggtitle("Branch Site Mutant Pilot vs Matched Pooled Oligo") +
  labs(x = "Nucleotide, Position",
       y = "Normalized SHAPE Reactivity")

### pearson vs One_144 == .547
struc <- shape %>%
  filter(source == "disrupted_disrupted") %>%
  ggplot(aes(x = as.factor(Nucleotide), y = Norm_profile)) +
  geom_col(fill = "#f1b6da") +
  geom_col(
    data = shape %>% filter(source == "One_144"),
    aes(x = as.factor(Nucleotide), y = -Norm_profile),
    fill = "#b8e186"
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_vline(xintercept = as.factor(188), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = as.factor(208), color = "black", linetype = "dashed") +
  ggtitle("Structure Mutant Pilot vs Matched Pooled Oligo") +
  labs(x = "Nucleotide, Position",
       y = "Normalized SHAPE Reactivity")
struc

pdf(file = "/Desktop/Dissertation/chapter_three/figures/combined.pdf",
    width = 12, height = 8, paper = "letter")
ggarrange(wt, bps, struc, ncol = 1, nrow = 3)
dev.off()


shape %>%
  filter(source == "disrupted_disrupted") %>%
  ggplot(aes(x = as.factor(Nucleotide), y = Norm_profile)) +
  geom_col(fill = "blue") +
  geom_col(
    data = shape %>% filter(source == "One_144"),
    aes(x = as.factor(Nucleotide), y = -Norm_profile),
    fill = "red"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, size = 4)) +
  geom_vline(xintercept = as.factor(188), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = as.factor(208), color = "black", linetype = "dashed") +
  ggtitle("Structure Mutant Pilot vs Matched Pooled Oligo") +
  labs(x = "Nucleotide, Position",
       y = "Normalized SHAPE Reactivity")



shape %>%
  filter(source %in% c("base_base", "bp_bp", "disrupted_disrupted",
                       "One_181", "One_144", "One_205")) %>%
  group_by(rep, source) %>%
  summarise(mean = mean(Norm_profile, na.rm = TRUE),
            median = median(Norm_profile, na.rm = T),
            sd = sd(Norm_profile, na.rm = TRUE))
