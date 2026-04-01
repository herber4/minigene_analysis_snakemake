library(tidyverse)
library(ggpubr)

shape <- read.table(file = "/Users/herber4/Desktop/Dissertation/chapter_three/data/NIA_Pilot_vs_Pooled_Shape_Data_For_Figures.txt",
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


pdf(file = "/Users/herber4/Desktop/Dissertation/chapter_three/figures/NIA_Pilot_vs_Pool/WT_Pilot_vs_WT_Pooled.pdf",
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

wt

wt_pool <- shape %>%
  filter(str_detect(source, "One_205|base_base"))

library(dplyr)
library(tidyr)

df_wide <- wt_pool %>%
  select(Nucleotide, Norm_profile, source) %>%
  pivot_wider(names_from = source, values_from = Norm_profile)

cor_test <- cor.test(
  df_wide$base_base,
  df_wide$One_205,
  method = "spearman",
  exact = FALSE
)

cor_test
rho <- cor_test$estimate
p_val <- cor_test$p.value

n <- sum(complete.cases(df_wide))

z <- atanh(rho)              # Fisher transform
se_z <- 1 / sqrt(n - 3)      # SE in z-space
se_rho <- (1 - rho^2) * se_z # delta method back-transform

z_obs <- atanh(rho)
z_null <- atanh(0.999999)  # can't use exactly 1, use ~1

z_score <- (z_obs - z_null) / se_z
p_val_vs1 <- 2 * pnorm(-abs(z_score))

list(
  spearman_rho = rho,
  p_value_vs_0 = p_val,
  standard_error = se_rho,
  p_value_vs_1 = p_val_vs1
)
