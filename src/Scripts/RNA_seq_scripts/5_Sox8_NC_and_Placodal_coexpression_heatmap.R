#   Script written by Jake Leese, 16th August 2024.
#   Purpose:
# - Plotting heatmap of NC_and_Placodal Sox8 coexpression_correlation 

setwd("/data/scRNAseq_objects")
.libPaths("/R/libs/R_scvi_integration_v3.4")

library(Seurat)
library(usethis)
library(devtools)
library(TFBSTools)
library(WGCNA)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

#### READING IN SOX8 CO-EXPRESSION VALUES AND FORMATTING FOR CLUSTERING / HEATMAP PLOTTING

# Read in .csv files. Read column 1 (genes) as row names, rather than a data column
Sox8_NC_TF_co_expression <- read.csv("Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages.csv", row.names = 1)

# Converting dataframe to numerical matrix in order to be read by complex heatmap         
correlation_matrix <- data.matrix(Sox8_NC_TF_co_expression)

# Identify rows whose sum values fall below a given threshold. Creates a logical vector which can then be used to filter out rows
Low_value_rows <- rowSums(abs(correlation_matrix)) > 0.5
correlation_matrix_filtered <- correlation_matrix[Low_value_rows, , drop = FALSE]

#### CUSTOM CLUSTERING FOR GROUPING GENES BASED ON PATTERNS OF CO-EXPRESSION ACROSS STAGES AND CELL TYPE

# hierarchical clustering using ward.D2 method. Distances calculated using euclidean method.
HCluster <- hclust(dist(correlation_matrix_filtered, method = "euclidean"), method = "ward.D2")

#### PLOTTING DATA USING COMPLEX HEATMAP

# Plotting complex heatmap
Heatmap <- Heatmap(correlation_matrix_filtered,
                   column_title = "Sox8 co-expression heatmap - NC and Placodal cells",
                   col = colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
                   heatmap_legend_param = list(title = "Correlation",
                                               at = c(-0.5, 0, 0.5),
                                               labels = c("-0.5", "0", "0.5")),
                   cluster_rows = HCluster,
                   cluster_columns = FALSE,
                   show_row_names = FALSE,
                   show_column_names = TRUE,
                   row_split = 9)

### EXTRACTING GENE CLUSTER IDs

# Show the heatmap
Heatmap <- draw(Heatmap)

# Extract clusters and check/confirm size of clusters. 
Cluster_list <- row_order(Heatmap)
lapply(Cluster_list, function(x) length(x))

# Loop to extract genes for each cluster
for (i in 1:length(Cluster_list)) {
  if (i == 1) {                                                                 # This if statement provides one iteration for the initial out object to be created without overwriting the updated version each time
    clu <- t(t(row.names(correlation_matrix_filtered[Cluster_list[[i]],]))) # extracts gene names that correspond to the row numbers of cluster, i. These row numbers are the same as the row numbers from the original matrix. The t(t()) functions transform these values (gene names) into a column ready to for cbind() 
    out <- cbind(clu, paste("Cluster", i, sep="_"))
    colnames(out) <- c("Gene", "Cluster")
  } else {
    clu <- t(t(row.names(correlation_matrix_filtered[Cluster_list[[i]],])))
    clu <- cbind(clu, paste("Cluster", i, sep="_"))
    out <- rbind(out, clu)
  }
}

# Save Gene clusters to excel file
write.csv(out, "Sox8_coexpression_heatmap_NC_and_Placodal_0.75_min_clusters.csv")
