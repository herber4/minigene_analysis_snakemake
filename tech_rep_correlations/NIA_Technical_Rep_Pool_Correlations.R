library(tidyverse)




rep_one <- "DownloadedSequenceData/austin/MAP_NovaSeq/nia_txsome_rep_one/shapemapper/A/profs"
rep_two <- "DownloadedSequenceData/austin/MAP_NovaSeq/nia_txsome_rep_two/shapemapper/A/profs"

rep_one <- list.files(rep_one, pattern = "_profile\\.txt$", full.names = TRUE)
rep_two <- list.files(rep_two, pattern = "_profile\\.txt$", full.names = TRUE)
# Get all *_profile.txt files


# Read and combine all files, adding a 'source' column with the filename
rep_one <- rep_one %>%
  map_df(~ read.table(.x, header = TRUE, sep = "\t") %>%
           mutate(source = basename(.x)))
rep_one$rep <- "One"
rep_two <- rep_two %>%
  map_df(~ read.table(.x, header = TRUE, sep = "\t") %>%
           mutate(source = basename(.x)))
rep_two$rep <- "Two"

rep_one$source <- sub("^((?:[^_]+_){1}[^_]+).*", "\\1", rep_one$source)
rep_two$source <- sub("^((?:[^_]+_){1}[^_]+).*", "\\1", rep_two$source)

one_filt <- rep_one %>%
  filter(Nucleotide >= 56,
         Nucleotide <= 225) %>%
  mutate(
    #Norm_profile = ifelse(is.nan(Norm_profile), 0, Norm_profile),
    Norm_profile = ifelse(Norm_profile < 0, 0, Norm_profile),
    Norm_profile = ifelse(Norm_profile > 4, 4, Norm_profile)
  )
two_filt <- rep_two %>%
  filter(Nucleotide >= 56,
         Nucleotide <= 225) %>%
  mutate(
    #Norm_profile = ifelse(is.nan(Norm_profile), 0, Norm_profile),
    Norm_profile = ifelse(Norm_profile < 0, 0, Norm_profile),
    Norm_profile = ifelse(Norm_profile > 4, 4, Norm_profile)
  )



one_wide <- one_filt[,c(1,28,30)]
one_wide <- one_wide %>%
  pivot_wider(names_from = source, values_from = Norm_profile, names_prefix = "Rep_One_")
two_wide <- two_filt[,c(1,28,30)]
two_wide <- two_wide %>%
  pivot_wider(names_from = source, values_from = Norm_profile, names_prefix = "Rep_Two_")

wide <- merge(one_wide, two_wide, by.x = "Nucleotide")

sample_cols <- setdiff(names(wide), "Nucleotide")

# --- 2) compute Spearman correlation matrix (pairwise, using available pairs) ---
mat <- wide %>%
  select(all_of(sample_cols)) %>%
  as.matrix()

spearman_mat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")

# --- 3) convert upper triangle of matrix to tidy long format ---
idx <- which(upper.tri(spearman_mat), arr.ind = TRUE)
spearman_pairs <- tibble(
  sample1 = rownames(spearman_mat)[idx[, "row"]],
  sample2 = colnames(spearman_mat)[idx[, "col"]],
  spearman = spearman_mat[idx]
) %>% arrange(sample1, sample2)

# inspect
spearman_pairs
pearsan_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")

# --- 3) convert upper triangle of matrix to tidy long format ---
idx <- which(upper.tri(pearsan_mat), arr.ind = TRUE)
pearsan_pairs <- tibble(
  sample1 = rownames(pearsan_mat)[idx[, "row"]],
  sample2 = colnames(pearsan_mat)[idx[, "col"]],
  pearson = pearsan_mat[idx]
) %>% arrange(sample1, sample2)

spearman_pairs$source <- "Spearman"
colnames(spearman_pairs)[3] <- "correlation"
pearsan_pairs$source <- "Pearson"
colnames(pearsan_pairs)[3] <- "correlation"
corrs <- rbind(pearsan_pairs, spearman_pairs)

ggplot(corrs, aes(y = correlation, x = source, group = source)) +
  geom_boxplot()

filt_corrs <- corrs %>% 
  subset(str_detect(sample1, "Rep_One") & str_detect(sample2, "Rep_Two"))


filt_corrs %>%
  filter(sample1 %in% c("Rep_One_One_205", "Rep_One_One_188",
                        "Rep_One_One_213", "Rep_One_One_73",
                        "Rep_One_One_26")) %>%
  ggplot(aes(x = sample1, y = correlation, fill = source)) +
  geom_boxplot()


filt_corrs <- filt_corrs %>%
  filter(str_extract(sample1, "[^_]+_[^_]+$") == str_extract(sample2, "[^_]+_[^_]+$"))

ggplot(filt_corrs, aes(y = correlation, x = source, group = source)) +
  geom_boxplot()

within_rep <- corrs %>% 
  subset(str_detect(sample1, "Rep_One_One") & str_detect(sample2, "Rep_Two_One"))

within_rep %>%
  filter(sample1 %in% c("Rep_One_One_205", "Rep_One_One_188",
                        "Rep_One_One_213", "Rep_One_One_73",
                        "Rep_One_One_26")) %>%
  ggplot(aes(x = sample1, y = correlation, fill = source)) +
  geom_boxplot()


btwn_reps <- corrs %>%
  subset(str_detect(sample1, "Rep_One_One") & str_detect(sample2, "Rep_Two_One"))

library(dplyr)

btwn_reps_subset <- btwn_reps %>%
  mutate(
    suffix1 = sub("^[^_]+_[^_]+_", "", sample1),
    suffix2 = sub("^[^_]+_[^_]+_", "", sample2)
  ) %>%
  filter(suffix1 == suffix2) %>%
  select(-suffix1, -suffix2)

pdf(file = "DownloadedSequenceData/austin/MAP_NovaSeq/R_scripts/NIA_Between_Reps_Corrlations.pdf",
    width = 12, height = 8, paper = "letter")
ggplot(btwn_reps_subset, aes(x = source, y = correlation)) +
  geom_boxplot() +
  theme_bw() +
  ylim(0, 1)
dev.off()
