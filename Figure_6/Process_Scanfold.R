library(ggplot2)
library(dplyr)
library(stringr)
library(dplyr)

file_list <- list.files("/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/structure_similarity/scanfold/base_pairs", pattern = "\\.txt$", full.names = TRUE)  # List of file names in the directory

combined_data <- data.frame()  # Initialize an empty dataframe to store the combined data

for (file in file_list) {
  data <- read.table(file, header = FALSE, skip = 1)  # Read the file
  gene_name <- tools::file_path_sans_ext(basename(file))  # Extracting the file name without extension
  data$gene <- gene_name  # Adding a column named 'gene' with the file name
  combined_data <- bind_rows(combined_data, data)  # Append data to the combined dataframe
}

combined_data <- rename(combined_data,
                        nt = V1,
                        i = V2,
                        j = V3,
                        avgMFE = V4,
                        avgZ = V5,
                        avgED = V6)


# Now `combined_data` contains all the data from the individual files with a 'gene' column
combined_data$nt <- gsub("[^0-9]", "", combined_data$nt)
combined_data$Nucleotide <- combined_data$nt

write.table(combined_data, file = "/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/structure_similarity/scanfold_per_base_Z_scores.txt",
            sep = "\t", row.names = F, quote = F)
