# Coexpression calculations and clustering analysis and heatmap generation of HH7 and HH8 stages only
# Based off "4_Calculating_Sox8_*_correlations" and "5_Sox8_coexpression_clustering_and_heatmap.R" scripts

setwd("/data//Sox8_binding_partner_analysis/scRNAseq_objects")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(TFBSTools)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)

source("/data/Sox8_binding_partner_analysis/Functions/Extract_HM_Clusters.R")


#### CALCULATE COEXPRESSION MATRICES (skip this step once it's been done once) ####

# Load and check the Seurat staged objects for NC and Placodal cells. Using those with Transcription factors only
HH7_NC_TFs <- readRDS("HH7_NC_TFs")
HH8_NC_TFs <- readRDS("HH8_NC_TFs")
HH7_Placodal_TFs <- readRDS("HH7_Placodal_TFs")
HH8_Placodal_TFs <- readRDS("HH8_Placodal_TFs")

# extract count data from NC_only and placode_only SEURAT object, and then subset based on TF genes only
HH7_NC_TFs_counts <- t(HH7_NC_TFs[["SCT"]]$counts)
HH8_NC_TFs_counts <- t(HH8_NC_TFs[["SCT"]]$counts)
HH7_Placodal_TFs_counts <- t(HH7_Placodal_TFs[["SCT"]]$counts)
HH8_Placodal_TFs_counts <- t(HH8_Placodal_TFs[["SCT"]]$counts)

# calculate correlation matrices
HH7_NC_TF_correlation_matrix <- WGCNA::cor(HH7_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH7_NC <- data.frame(as.list(HH7_NC_TF_correlation_matrix["SOX8", ]))

HH8_NC_TF_correlation_matrix <- WGCNA::cor(HH8_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH8_NC <- data.frame(as.list(HH8_NC_TF_correlation_matrix["SOX8", ]))

HH7_Placodal_TF_correlation_matrix <- WGCNA::cor(HH7_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH7_Placodal <- data.frame(as.list(HH7_Placodal_TF_correlation_matrix["SOX8", ]))

HH8_Placodal_TF_correlation_matrix <- WGCNA::cor(HH8_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH8_Placodal <- data.frame(as.list(HH8_Placodal_TF_correlation_matrix["SOX8", ]))

# Make a dataframe from Sox8 Spearman's rank values for each stage
# First, assign rownames to each Sox8 correlation dataframe.
rownames(Sox8_co_expressed_TFs_HH7_NC) <- "HH7"
rownames(Sox8_co_expressed_TFs_HH8_NC) <- "HH8"

rownames(Sox8_co_expressed_TFs_HH7_Placodal) <- "HH7"
rownames(Sox8_co_expressed_TFs_HH8_Placodal) <- "HH8"

# Then combine the dataframes into 1 using rbind()
Sox8_NC_co_expressed_TFs_HH7_HH8 <- rbind(Sox8_co_expressed_TFs_HH7_NC,
                                          Sox8_co_expressed_TFs_HH8_NC)
Sox8_NC_co_expressed_TFs_HH7_HH8 <- t(Sox8_NC_co_expressed_TFs_HH7_HH8)

Sox8_Placodal_co_expressed_TFs_HH7_HH8 <- rbind(Sox8_co_expressed_TFs_HH7_Placodal,
                                                Sox8_co_expressed_TFs_HH8_Placodal)
Sox8_Placodal_co_expressed_TFs_HH7_HH8 <- t(Sox8_Placodal_co_expressed_TFs_HH7_HH8)

# Remove SOX8 from dataframe and Replace all NA values with a 0 and save to disk
Sox8_NC_co_expressed_TFs_HH7_HH8 <- subset(Sox8_NC_co_expressed_TFs_HH7_HH8, rownames(Sox8_NC_co_expressed_TFs_HH7_HH8) != "SOX8")
Sox8_NC_co_expressed_TFs_HH7_HH8[is.na(Sox8_NC_co_expressed_TFs_HH7_HH8)] <- 0

write.csv(Sox8_NC_co_expressed_TFs_HH7_HH8, "CoExpression_values/Sox8_NC_co_expressed_TFs_HH7_HH8.csv")

Sox8_Placodal_co_expressed_TFs_HH7_HH8 <- subset(Sox8_Placodal_co_expressed_TFs_HH7_HH8, rownames(Sox8_Placodal_co_expressed_TFs_HH7_HH8) != "SOX8")
Sox8_Placodal_co_expressed_TFs_HH7_HH8[is.na(Sox8_Placodal_co_expressed_TFs_HH7_HH8)] <- 0

write.csv(Sox8_Placodal_co_expressed_TFs_HH7_HH8, "CoExpression_values/Sox8_Placodal_co_expressed_TFs_HH7_HH8.csv")


#### READING IN SOX8 CO-EXPRESSION VALUES AND FORMATTING FOR CLUSTERING / HEATMAP PLOTTING

# Read in .csv files. Read column 1 (genes) as row names, rather than a data column
Sox8_NC_TF_HH7_HH8_co_expression <- read.csv("CoExpression_values/Sox8_NC_co_expressed_TFs_HH7_HH8.csv", row.names = 1)
Sox8_Placodal_TF_HH7_HH8_co_expression <- read.csv("CoExpression_values/Sox8_Placodal_co_expressed_TFs_HH7_HH8.csv", row.names = 1)

# Rename columns in dataframes to represent stage_celltype
colnames(Sox8_NC_TF_HH7_HH8_co_expression) <- paste(colnames(Sox8_NC_TF_HH7_HH8_co_expression), "NC", sep = "_")
colnames(Sox8_Placodal_TF_HH7_HH8_co_expression) <- paste(colnames(Sox8_Placodal_TF_HH7_HH8_co_expression), "Placodal", sep = "_")

# Combine dataframes and reorder columns 
correlation_matrix_HH7_HH8 <- cbind(Sox8_NC_TF_HH7_HH8_co_expression, Sox8_Placodal_TF_HH7_HH8_co_expression)
correlation_matrix_HH7_HH8 <- correlation_matrix_HH7_HH8[, c(1,3,2,4)] # re-ordering columns

# Filtering out all of the rows in the matrix that do not meet a minimum correlation/anti-correlation threshold
# Identify rows whose sum values fall below a given threshold. Creates a logical vector which can then be used to filter out rows
Low_value_rows <- rowSums(abs(correlation_matrix_HH7_HH8)) > 0.2
correlation_matrix_HH7_HH8_filtered <- correlation_matrix_HH7_HH8[Low_value_rows, , drop = FALSE]


#### PLOTTING DATA USING COMPLEX HEATMAP

# hierarchical clustering using ward.D2 method. Distances calculated using euclidean method.
HCluster <- hclust(dist(correlation_matrix_HH7_HH8_filtered, method = "euclidean"), method = "ward.D2")

# Converting dataframe to numerical matrix in order to be read by complex heatmap         
correlation_matrix_filtered_mat <- data.matrix(correlation_matrix_HH7_HH8_filtered)

# defining a column split vector. unique(col_split) then respects the ordering of the input vector (https://stackoverflow.com/questions/55799509/avoid-re-ordering-of-rows-and-columns-in-a-heatmap-r)
col_split <- rep(c("HH7", "HH8"), each = 2)

# Creating heatmap column labels and defining colours
stages <- c("HH7", "HH7", "HH8", "HH8")
CellType <- c("NC", "Placodal", "NC", "Placodal")

annotation_col <- data.frame(
  Stage = stages,
  CellType = CellType)
rownames(annotation_col) <- colnames(correlation_matrix_HH7_HH8_filtered)

# Define colors for annotation
annotation_colors <- list(
  Stage = c(HH7 = "lightblue", HH8 = "cyan4"),
  CellType = c(NC = "pink", Placodal = "aquamarine3")
)

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
                   row_split = 20,
                   column_split = factor(col_split, levels = unique(col_split))
)

# Show the heatmap
Heatmap <- draw(Heatmap)

#### EXTRACTING GENE CLUSTER IDs ####

# Custom function to extract cluster gene IDs from a heatmap
Sox8_coexpression_clusters_HH7_HH8_min_0.2_20C <- Extract_HM_Clusters(Heatmap, correlation_matrix_filtered_mat)

# Save Gene clusters to excel file
write.csv(Sox8_coexpression_clusters_HH7_HH8_min_0.2_20C, "Heatmap_clusters/Sox8_coexpression_clusters_HH7_HH8_min_0.2_20C.csv")
