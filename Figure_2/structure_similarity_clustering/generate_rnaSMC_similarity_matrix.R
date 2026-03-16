library(ggplot2)
library(RNAsmc)

setwd("/DownloadedSequenceData/USER/MAP_NovaSeq/structure_similarity/reps_combined")

# List all files in the directory with a .ct extension
ct_files <- list.files(pattern = "\\.ct$")

# Create an empty list to store data frames
ct_dataframes <- list()

# Loop through the .ct files and read them into data frames
for (file in ct_files) {
  # Read the file into a data frame (adjust read.csv parameters as needed)
  df <- read.table(file, skip = 1)
  
  # Store the data frame in the list with the same name as the file
  ct_dataframes[[file]] <- df
}


clus <- RNAstrCluster(ct_dataframes)

tmp <- as.data.frame(clus[["simility_mat"]])
tmp$samps <- rownames(tmp)
rownames(tmp) <- NULL

write.table(tmp, file = "/DownloadedSequenceData/USER/MAP_NovaSeq/structure_similarity/combined_reps_similarity_matrix.txt",
            sep = "\t", quote = F)
