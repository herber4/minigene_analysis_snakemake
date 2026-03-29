library(tidyverse)
library(stringr)
setwd("/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts/")

# all files
files <- c(
  "rep_one_k700e_samtools.depth.counts.txt",
  "rep_one_WT_samtools.depth.counts.txt",
  "rep_three_k700e_samtools.depth.counts.txt",
  "rep_three_WT_samtools.depth.counts.txt",
  "rep_two_k700e_samtools.depth.counts.txt",
  "rep_two_WT_samtools.depth.counts.txt"
)

# helper function to process one file
process_file <- function(file_path) {
  df <- read.table(file_path, sep = "\t", header = TRUE)
  
  # remove BAM suffix
  colnames(df) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(df))
  
  # pivot longer
  df_long <- df %>%
    pivot_longer(
      cols = -c(CHROM, POS),
      names_to = "sample",
      values_to = "coverage"
    ) %>%
    group_by(POS, sample) %>%
    summarise(coverage = sum(coverage), .groups = "drop")
  
  # create column name based on file
  # extract replicate and K700E/WT info
  file_base <- basename(file_path)
  file_base <- gsub("_samtools.depth.counts.txt", "", file_base)
  
  # convert to MT/WT naming
  if (grepl("k700e", file_base, ignore.case = TRUE)) {
    col_name <- paste0("MT_", str_to_title(gsub("rep_", "Rep_", str_extract(file_base, "rep_[a-z]+"))))
  } else {
    col_name <- paste0("WT_", str_to_title(gsub("rep_", "Rep_", str_extract(file_base, "rep_[a-z]+"))))
  }
  
  df_long <- df_long %>%
    rename(!!col_name := coverage)
  
  return(df_long)
}

# apply to all files
list_dfs <- lapply(files, process_file)

# merge by POS and sample
final_df <- Reduce(function(x, y) full_join(x, y, by = c("POS", "sample")), list_dfs)
rm(list_dfs, mt_long, mt_valid, tmp, mt_collapsed, mt)

feats <- read.table(file = "annotation_features.txt",
                    sep = "\t", header = TRUE)

mas <- merge(feats, final_df, by = "sample")


# Identify the columns with counts (everything except POS and sample)
count_cols <- setdiff(colnames(final_df), c("POS", "sample"))

# Compute CPM
cpm_df <- final_df %>%
  mutate(across(all_of(count_cols), 
                ~ . / sum(.) * 1e6,  # counts divided by column sum * 1 million
                .names = "{.col}_CPM"))

#write.table(cpm_df, file = "cpm_df.txt",
#            sep = "\t", quote = F)
cpm_df <- cpm_df[,c(1:2, 9:14)]

cpm_long <- cpm_df %>%
  pivot_longer(cols = 3:8,
               values_to = "cpm", names_to = "sample_name")

cpm_long <- merge(feats, cpm_long, by = "sample")

cpm_long %>%
  filter(POS > 359, POS < 383,
         annot == "NAGNAG") %>%
  ggplot(aes(x = factor(POS), y = cpm, fill = sample_name)) +
  geom_boxplot()
  
#### Geom Tile Plot 
{
library(ggplot2)
library(dplyr)

cpm_long %>%
  filter(POS > 359, POS < 385) %>%
  ggplot(aes(x = POS, y = sample, fill = cpm)) +
  geom_tile() +
  scale_fill_viridis_c(name = "CPM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Position",
    y = "Sample",
    title = "CPM Heatmap (POS 360–382)"
  ) +
  facet_wrap(~sample_name)
}
#####

#c3ss vs 3ss

library(dplyr)

ss_usage <- cpm_long %>%
  filter(POS %in% c(362, 363, 382, 383)) %>%
  mutate(
    site = case_when(
      POS %in% c(362, 363) ~ "C3SS",
      POS %in% c(382, 383) ~ "3SS"
    )
  ) %>%
  group_by(sample, sample_name, site) %>%
  summarise(
    site_cpm = sum(cpm, na.rm = TRUE),
    .groups = "drop"
  )

ss_ratio <- ss_usage %>%
  tidyr::pivot_wider(
    names_from = site,
    values_from = site_cpm
  ) %>%
  mutate(
    C3SS_ratio = C3SS / (C3SS + `3SS`)
  )


ss_ratio <- ss_ratio %>%
  mutate(
    genotype = if_else(grepl("^MT", sample_name), "MT", "WT")
  )

ss_summary <- ss_ratio %>%
  group_by(sample, genotype) %>%
  summarise(
    mean_ratio = mean(C3SS_ratio, na.rm = TRUE),
    sd_ratio   = sd(C3SS_ratio, na.rm = TRUE),
    .groups = "drop"
  )

ss_summary$sd_ratio <- NULL 
ss_compare <- ss_summary %>%
  pivot_wider(names_from = genotype, values_from = mean_ratio)


stats <- ss_ratio %>%
  group_by(sample) %>%
  summarise(
    p_value = t.test(
      C3SS_ratio ~ genotype
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH")
  )

final_results <- ss_compare %>%
  left_join(stats, by = "sample")


ggplot(ss_ratio, aes(genotype, C3SS_ratio, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  facet_wrap(~ sample, scales = "free_y") +
  theme_minimal() +
  labs(y = "C3SS / (C3SS + 3SS)")
