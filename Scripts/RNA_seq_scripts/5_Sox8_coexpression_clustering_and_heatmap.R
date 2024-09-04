#   Script written by Jake Leese, 14th August 2024.
#   Purpose:
# - Combining NC and Placodal Sox8 correlation coefficients dataframes into a single matrix  
# - Hierarchhical clustering of genes based on patterns of Sox8 co-expression across stages and tissue type
# - Plotting everything on a heatmap, with genes ordered into clusters, or "modules"

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
Sox8_NC_TF_co_expression <- read.csv("Sox8_NC_co_expressed_TFs_combined_stages.csv", row.names = 1)
Sox8_Placodal_TF_co_expression <- read.csv("Sox8_Placodal_co_expressed_TFs_combined_stages.csv", row.names = 1)

# Rename columns in dataframes to represent stage_celltype
colnames(Sox8_NC_TF_co_expression) <- paste(colnames(Sox8_NC_TF_co_expression), "NC", sep = "_")
colnames(Sox8_Placodal_TF_co_expression) <- paste(colnames(Sox8_Placodal_TF_co_expression), "Placodal", sep = "_")

# Combine dataframes and reorder columns 
correlation_matrix <- cbind(Sox8_NC_TF_co_expression, Sox8_Placodal_TF_co_expression)
correlation_matrix <- correlation_matrix[, c(1,7,2,8,3,9,4,10,5,11,6,12)] # re-ordering columns

# Filtering out all of the rows in the matrix that do not meet a minimum correlation/anti-correlation threshold
# Identify rows whose sum values fall below a given threshold. Creates a logical vector which can then be used to filter out rows
Low_value_rows <- rowSums(abs(correlation_matrix)) > 0.75
correlation_matrix_filtered <- correlation_matrix[Low_value_rows, , drop = FALSE]

#### CUSTOM CLUSTERING FOR GROUPING GENES BASED ON PATTERNS OF CO-EXPRESSION ACROSS STAGES AND CELL TYPE

# hierarchical clustering using ward.D2 method. Distances calculated using euclidean method.
HCluster <- hclust(dist(correlation_matrix_filtered, method = "euclidean"), method = "ward.D2")


#### CREATING ANNOTATION AND ANNOTATION COLOUR OBJECTS FOR HEATMAPS

# Creating heatmap column labels
stages <- c("HH7", "HH7", "HH8", "HH8", "HH9", "HH9", "HH12", "HH12", "HH14", "HH14", "HH16", "HH16")
CellType <- c("NC", "Placodal", "NC", "Placodal", "NC", "Placodal", "NC", "Placodal", "NC", "Placodal", "NC", "Placodal")

annotation_col <- data.frame(
  Stage = stages,
  CellType = CellType)
rownames(annotation_col) <- colnames(correlation_matrix_filtered)

# Define colors for annotation
annotation_colors <- list(
  Stage = c(HH7 = "lightblue", HH8 = "cyan4", HH9 = "palegreen",  HH12 = "gold1", HH14 = "lightsalmon", HH16 = "tomato"),
  CellType = c(NC = "pink", Placodal = "aquamarine3")
)

#### PLOTTING HEATMAP USING PHEATMAP

Sox8_co_expression_heatmap <- 
  pheatmap::pheatmap(correlation_matrix_filtered,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = FALSE,
         annotation_legend = TRUE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         cutree_rows = 10,
         color = colorRampPalette(c("blue", "white", "red"))(25),
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         gaps_col = c(2,4,6,8,10),
         display_numbers = FALSE,
         fontsize_row = 2,
         main = "Sox8 Co-expression Heatmap")
        
# set white color for dendrogram lines and text

#### PLOTTING DATA USING COMPLEX HEATMAP

# Converting dataframe to numerical matrix in order to be read by complex heatmap         
correlation_matrix_filtered_mat <- data.matrix(correlation_matrix_filtered)

# defining a column split vector. unique(col_split) then respects the ordering of the input vector (https://stackoverflow.com/questions/55799509/avoid-re-ordering-of-rows-and-columns-in-a-heatmap-r)
col_split <- rep(c("HH7", "HH8", "HH9", "HH12", "HH14", "HH16"), each = 2)

# Plotting complex heatmap
Heatmap <- Heatmap(correlation_matrix_filtered_mat,
        top_annotation = HeatmapAnnotation(
          Stage = annotation_col$Stage,
          col = list(Stage = annotation_colors$Stage),
          show_legend = FALSE,
          show_annotation_name = FALSE),
        bottom_annotation = HeatmapAnnotation(
          CellType = annotation_col$CellType,
          col = list(CellType = annotation_colors$CellType),
          show_annotation_name = FALSE),
        col = colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
        heatmap_legend_param = list(title = "Correlation",
                                    at = c(-0.5, 0, 0.5),
                                    labels = c("-0.5", "0", "0.5")),
        cluster_rows = HCluster,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        row_split = 8,
        column_split = factor(col_split, levels = unique(col_split))
        )

### EXTRACTING GENE CLUSTER IDs

# Show the heatmap
Heatmap <- draw(Heatmap)

# Extract clusters and check/confirm size of clusters. 
Cluster_list <- row_order(Heatmap)
lapply(Cluster_list, function(x) length(x))

# Loop to extract genes for each cluster
for (i in 1:length(Cluster_list)) {
  if (i == 1) {                                                                 # This if statement provides one iteration for the initial out object to be created without overwriting the updated version each time
    clu <- t(t(row.names(correlation_matrix_filtered_mat[Cluster_list[[i]],]))) # extracts gene names that correspond to the row numbers of cluster, i. These row numbers are the same as the row numbers from the original matrix. The t(t()) functions transform these values (gene names) into a column ready to for cbind() 
    out <- cbind(clu, paste("Cluster", i, sep="_"))
    colnames(out) <- c("Gene", "Cluster")
  } else {
    clu <- t(t(row.names(correlation_matrix_filtered_mat[Cluster_list[[i]],])))
    clu <- cbind(clu, paste("Cluster", i, sep="_"))
    out <- rbind(out, clu)
  }
}

# Save Gene clusters to excel file
write.csv(out, "Sox8_co_expression_clusters_0.75_min.csv")
