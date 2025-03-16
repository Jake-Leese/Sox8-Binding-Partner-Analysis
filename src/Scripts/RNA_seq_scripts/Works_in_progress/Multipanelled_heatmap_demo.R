# Example matrix (replace this with your actual data)
data_matrix <- matrix(runif(50), nrow=5)
rownames(data_matrix) <- paste0("Gene", 1:5)
colnames(data_matrix) <- c("Stage1_Sample1", "Stage1_Sample2", "Stage2_Sample1", "Stage2_Sample2", "Stage3_Sample1", "Stage3_Sample2",
                           "Stage4_Sample1", "Stage4_Sample2", "Stage5_Sample1", "Stage5_Sample2")

# Create a vector for stages and samples
stages <- c("Stage1", "Stage1", "Stage2", "Stage2", "Stage3", "Stage3", "Stage4", "Stage4", "Stage5", "Stage5")
samples <- c("Sample1", "Sample2", "Sample1", "Sample2", "Sample1", "Sample2", "Sample1", "Sample2", "Sample1", "Sample2")

# Create a data frame for column annotation
annotation_col <- data.frame(
  Stage = stages,
  Sample = samples
)
rownames(annotation_col) <- colnames(data_matrix)

# Define colors for annotation
annotation_colors <- list(
  Stage = c(Stage1 = "lightblue", Stage2 = "lightgreen", Stage3 = "lightpink", Stage4 = "lightyellow1", Stage5 = "lightsalmon"),
  Sample = c(Sample1 = "purple", Sample2 = "orange")
)

# Plot the heatmap
pheatmap(data_matrix,
         annotation_col = annotation_col,  # Add column annotations
         annotation_colors = annotation_colors,  # Add annotation colors (optional)
         cluster_rows = TRUE,  # Cluster genes (rows)
         cluster_cols = FALSE,  # Disable clustering for columns to maintain order
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show sample names
         fontsize = 10,  # Adjust text size
         main = "Gene Expression Heatmap")

plot(heatmap)
