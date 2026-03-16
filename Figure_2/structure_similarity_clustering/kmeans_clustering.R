library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(factoextra)
setwd("")
df <- read.table(file = "combined_reps_similarity_matrix.txt",
                 sep = "\t", header = TRUE)


df <- df[,c(499, 1:498)]



mat <- df
colnames(mat) <- gsub(".ct", "", colnames(mat))
rownames(mat) <- mat$samps
mat$samps <- NULL
rownames(mat) <- gsub(".ct", "", rownames(mat))

mat <- as.matrix(mat)

## ------------------------------------------------------------------
## Convert rnaSMC similarity matrix to a distance matrix
## ------------------------------------------------------------------

# Convert input object to a numeric matrix
# rnaSMC outputs pairwise similarity scores ranging from 0–10,
# where 10 indicates identical RNA structures
sim <- as.matrix(mat)

# Convert similarity to distance so that:
#   high similarity → small distance
#   low similarity  → large distance
# This linear inversion preserves relative ordering
dist_mat <- 10 - sim


## ------------------------------------------------------------------
## Optional normalization of distances
## ------------------------------------------------------------------

# Scale distances to the range [0, 1]
# This is not strictly required but helps stabilize downstream methods
# and makes distances comparable across datasets
dist_mat <- dist_mat / max(dist_mat)

# Ensure that self-distances are exactly zero
# (important for distance-based methods)
diag(dist_mat) <- 0


## ------------------------------------------------------------------
## Multidimensional scaling (MDS)
## ------------------------------------------------------------------

# Classical MDS embeds the distance matrix into Euclidean space
# k = 2 produces two coordinates per RNA structure for visualization
mds <- cmdscale(as.dist(dist_mat), k = 2)

# Convert MDS output to a data frame for compatibility with factoextra
mds_df <- as.data.frame(mds)

# Assign intuitive axis names
colnames(mds_df) <- c("Dim1", "Dim2")


## ------------------------------------------------------------------
## Determine an appropriate number of clusters (k)
## ------------------------------------------------------------------

# Elbow method: plots within-cluster sum of squares vs. k
# Look for a point where improvement begins to plateau
fviz_nbclust(mds_df, kmeans, method = "wss") +
  theme_minimal()

# Silhouette method: measures cluster separation quality
# Higher average silhouette width indicates better-defined clusters
fviz_nbclust(mds_df, kmeans, method = "silhouette") +
  theme_minimal()


## ------------------------------------------------------------------
## Run k-means clustering
## ------------------------------------------------------------------

# Set random seed for reproducibility
set.seed(123)

# Specify number of clusters (chosen from biological insight
# and/or the diagnostic plots above)
k <- 3   # example value

# Perform k-means clustering on the MDS coordinates
# nstart = 50 runs the algorithm multiple times to avoid local minima
km <- kmeans(mds_df, centers = k, nstart = 50)


## ------------------------------------------------------------------
## Save k-means clustering visualization to PDF
## ------------------------------------------------------------------

# Open PDF graphics device
pdf(file = "kmeans_clustering.pdf",
    width = 10, height = 10, paper = "letter")

# Visualize clusters using factoextra
# Points = RNA structures
# Ellipses = multivariate normal approximation of clusters
# Cluster centers are shown
fviz_cluster(
  km,
  data = mds_df,
  geom = "point",
  ellipse.type = "norm",
  repel = TRUE,
  show.clust.cent = TRUE
) +
  labs(
    title = "k-means clustering of RNA structures (rnaSMC)",
    x = "MDS dimension 1",
    y = "MDS dimension 2"
  ) +
  theme_minimal()

# Close PDF device and write file to disk
dev.off()


## ------------------------------------------------------------------
## Alternative on-screen visualization with different styling
## ------------------------------------------------------------------

# Produce the same clustering plot with a different color palette
# and ggplot theme (useful for exploratory work)
fviz_cluster(
  km,
  data = mds_df,
  palette = "jco",      # alternative palettes: "npg", "lancet"
  ggtheme = theme_bw(),
  pointsize = 2,
  alpha = 0.8
)


## ------------------------------------------------------------------
## Extract cluster assignments for downstream analysis
## ------------------------------------------------------------------

# Create a tidy data frame linking each RNA to its assigned cluster
# Useful for merging with metadata or structure annotations
cluster_assignments <- data.frame(
  RNA = rownames(mat),
  cluster = km$cluster
)


## ------------------------------------------------------------------
## Assess distance preservation from MDS
## ------------------------------------------------------------------

# Extract eigenvalues from classical MDS
# Larger positive eigenvalues indicate dimensions that
# preserve more of the original distance information
cmdscale(as.dist(dist_mat), eig = TRUE)$eig
