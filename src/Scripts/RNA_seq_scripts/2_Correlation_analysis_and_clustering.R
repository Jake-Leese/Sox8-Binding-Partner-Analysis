setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(TFBSTools)
library(JASPAR2024)
library(dplyr)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(universalmotif)
library(ggplot2)
library(ggrepel)
library(tidyr)

source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_HM_Clusters.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/subset_seurat_features.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_count_data.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Correlation_analysis.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/plot_correlation_against_counts.R")

# Load the original Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

# Subset Seurat object using custom subset_seurat_features() function that takes arguments
# Seurat_obj, DB ("JASPAR2024" (vertebrates) or "AnimalTFDB" (Ggal)) 

HHall_ectoderm_TFs <- subset_seurat_features(HHall_ectoderm, DB = "JASPAR2024")

#### EXTRACT COUNT DATA ####

# Custom function to extract count data for: Stage, Cell_type (ectoderm type), Sample (origin.ident)(optional), and Assay ("RNA" or "SCT"), data_type ("counts" or "data" (log.normalised)), and subset_cells_by_gene ("gene.name") 
HH8_placode_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH8", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
HH8_NC_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH8", Cell_type = "NC", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
HH9_placode_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
HH9_NC_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "NC", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")


#### CALCULATE CORRELATION COEFFICIENTS ####

# Custom function to calculate correlation coefficients 
# Function takes arguments: remove_null_count_genes (logical TRUE or FALSE), Test ("spearman" or "pearson"), Gene (Specify a gene to take all values against. If left blank, full correlation matrix will be returned)

# Returns correlation values for SOX8 only
HH8_placode_RNA_spearman_SOX8 <- Correlation_analysis(HH8_placode_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
HH8_NC_RNA_spearman_SOX8 <- Correlation_analysis(HH8_NC_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
HH9_placode_RNA_spearman_SOX8 <- Correlation_analysis(HH9_placode_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
HH9_NC_RNA_spearman_SOX8 <- Correlation_analysis(HH9_NC_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")

# Combine correlation values to a single dataframe
HH8_placode_RNA_spearman_SOX8 <- tibble::rownames_to_column(as.data.frame(t(HH8_placode_RNA_spearman_SOX8)), "gene") %>% rename("HH8_Placode" = 2)
HH8_NC_RNA_spearman_SOX8 <- tibble::rownames_to_column(as.data.frame(t(HH8_NC_RNA_spearman_SOX8)), "gene") %>% rename("HH8_NC" = 2)
HH9_placode_RNA_spearman_SOX8 <- tibble::rownames_to_column(as.data.frame(t(HH9_placode_RNA_spearman_SOX8)), "gene") %>% rename("HH9_Placode" = 2)
HH9_NC_RNA_spearman_SOX8 <- tibble::rownames_to_column(as.data.frame(t(HH9_NC_RNA_spearman_SOX8)), "gene") %>% rename("HH9_NC" = 2)

# Combine all dataframes by "gene"
combined_df <- full_join(HH8_placode_RNA_spearman_SOX8, HH8_NC_RNA_spearman_SOX8, by = "gene") %>%
  full_join(HH9_placode_RNA_spearman_SOX8, by = "gene") %>%
  full_join(HH9_NC_RNA_spearman_SOX8, by = "gene") 
combined_df[is.na(combined_df)] <- 0
rownames(combined_df) <- combined_df[[1]] # Change the gene column into the rownames
combined_df <- combined_df[, -1]

# remove SOX8 from correlation matrix
combined_df <- subset(combined_df, rownames(combined_df) != "SOX8")

# Filtering out all of the rows (genes) in the matrix that do not have at least one Sox8 Co-expression value of 0.1 or greater
correlation_matrix_HH8_HH9_0.1 <- combined_df %>% filter(if_any(everything(), ~ . >=0.1))


#### PLOTTING DATA USING COMPLEX HEATMAP

# hierarchical clustering using ward.D2 method. Distances calculated using euclidean method.
HCluster <- hclust(dist(correlation_matrix_HH8_HH9_0.1, method = "euclidean"), method = "ward.D2")

# Converting dataframe to numerical matrix in order to be read by complex heatmap         
correlation_matrix_filtered_mat <- data.matrix(correlation_matrix_HH8_HH9_0.1)

# defining a column split vector. unique(col_split) then respects the ordering of the input vector (https://stackoverflow.com/questions/55799509/avoid-re-ordering-of-rows-and-columns-in-a-heatmap-r)
col_split <- rep(c("HH8", "HH9"), each = 2)

# Creating heatmap column labels and defining colours
stages <- c("HH8", "HH8", "HH9", "HH9")
CellType <- c("Placode", "NC", "Placode", "NC")

annotation_col <- data.frame(
  Stage = stages,
  CellType = CellType)
rownames(annotation_col) <- colnames(correlation_matrix_HH8_HH9_0.1)

# Define colors for annotation
annotation_colors <- list(
  Stage = c(HH8 = "lightblue", HH9 = "darkolivegreen2"),
  CellType = c(NC = "sienna3", Placode = "darkcyan")
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
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   row_split = 12,
                   column_split = factor(col_split, levels = unique(col_split))
)

# Show the heatmap
Heatmap <- draw(Heatmap)

svg("HH8_HH9_SOX8_Coexpression_Heatmap_RNA_assay.svg", width = 8, height = 6)
draw(Heatmap)
dev.off()

# Custom function to extract cluster gene IDs from a heatmap
Sox8_coexpression_clusters_HH8_HH9_0.1_12C <- Extract_HM_Clusters(Heatmap, correlation_matrix_filtered_mat)

# Adding a cluster label column to the original correlation_matrix
Sox8_coexpression_clusters_HH8_HH9_0.1_12C <- as.data.frame(Sox8_coexpression_clusters_HH8_HH9_0.1_12C[match(rownames(correlation_matrix_filtered_mat),Sox8_coexpression_clusters_HH8_HH9_0.1_12C[,1]), ])
correlation_matrix_filtered_df <- as.data.frame(correlation_matrix_filtered_mat)
correlation_matrix_filtered_df$Heatmap_Cluster <- c(Sox8_coexpression_clusters_HH8_HH9_0.1_12C[,2])

# re-ordering based on cluster number
correlation_matrix_filtered_df$Heatmap_Cluster <- as.numeric(sub(".*_(\\d+)$", "\\1", correlation_matrix_filtered_df$Heatmap_Cluster)) # Changing the cluster column to contain only numbers. This allows us to re-order the dataframe based on the order of the clusters
correlation_matrix_filtered_df <- correlation_matrix_filtered_df[order(correlation_matrix_filtered_df$Heatmap_Cluster),]

# Saving co-expression values together along with the cluster number and cell type annotations across stages
# write.csv(correlation_matrix_filtered_df, "Heatmap_clusters/Sox8_coexpression_clusters_HH8_HH9_JASPAR_RNA_assay_0.1_12C_SOX8_cells.csv")


#### FILTERING TFs ####

# Using the mean + standard deviation as a threshold for co-expression. Using pre-filtered values for calculations (no 0 values and no SOX8)
HH8_Placode_TFs <- rownames(combined_df[combined_df$HH8_Placode > (sd(as.numeric(HH8_placode_RNA_spearman_SOX8[HH8_placode_RNA_spearman_SOX8 < 0.99])) + mean(as.numeric(HH8_placode_RNA_spearman_SOX8[HH8_placode_RNA_spearman_SOX8 < 0.99]))),]) # threshold = 0.109216
HH8_NC_TFs <- rownames(combined_df[combined_df$HH8_NC > (sd(as.numeric(HH8_NC_RNA_spearman_SOX8[HH8_NC_RNA_spearman_SOX8 < 0.99])) + mean(as.numeric(HH8_NC_RNA_spearman_SOX8[HH8_NC_RNA_spearman_SOX8 < 0.99]))),]) # threshold = 0.1171457
HH9_Placode_TFs <- rownames(combined_df[combined_df$HH9_Placode > (sd(as.numeric(HH9_placode_RNA_spearman_SOX8[HH9_placode_RNA_spearman_SOX8 < 0.99])) + mean(as.numeric(HH9_placode_RNA_spearman_SOX8[HH9_placode_RNA_spearman_SOX8 < 0.99]))),]) # threshold = 0.1486482
HH9_NC_TFs <- rownames(combined_df[combined_df$HH9_NC > (sd(as.numeric(HH9_NC_RNA_spearman_SOX8[HH9_NC_RNA_spearman_SOX8 < 0.99])) + mean(as.numeric(HH9_NC_RNA_spearman_SOX8[HH9_NC_RNA_spearman_SOX8 < 0.99]))),]) # threshold = 0.2083915 

# Pull out TFs that are over the threshold in both NC and Placode cells
HH8_NC_and_Placode_TFs <- intersect(HH8_Placode_TFs, HH8_NC_TFs)
HH9_NC_and_Placode_TFs <- intersect(HH9_Placode_TFs, HH9_NC_TFs)
HH8_Placode_TFs <- setdiff(HH8_Placode_TFs, HH8_NC_and_Placode_TFs)  # Remove TFs that are also highly correlated in NC
HH8_NC_TFs <- setdiff(HH8_NC_TFs, HH8_NC_and_Placode_TFs) # Remove TFs that are also highly correlated in Placode
HH9_Placode_TFs <- setdiff(HH9_Placode_TFs, HH9_NC_and_Placode_TFs)  
HH9_NC_TFs <- setdiff(HH9_NC_TFs, HH9_NC_and_Placode_TFs)

#### MAKING A HEATMAP USING JUST THOSE TFs SELECTED FOR DOWNSTREAM ANALYSIS

All_TFs <- unique(c(HH8_Placode_TFs, HH8_NC_TFs, HH8_NC_and_Placode_TFs, HH9_Placode_TFs, HH9_NC_TFs, HH9_NC_and_Placode_TFs)) # Collect all of the TFs from our vectors

# Subset combined_df based on All_TFs vector
correlation_matrix_HH8_HH9_Filtered <- combined_df[rownames(combined_df) %in% All_TFs, ]

# hierarchical clustering using ward.D2 method. Distances calculated using euclidean method.
HCluster <- hclust(dist(correlation_matrix_HH8_HH9_Filtered, method = "euclidean"), method = "ward.D2")

# Converting dataframe to numerical matrix in order to be read by complex heatmap         
correlation_matrix_filtered_mat <- data.matrix(correlation_matrix_HH8_HH9_Filtered)

rownames(annotation_col) <- colnames(correlation_matrix_HH8_HH9_Filtered)

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
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   row_split = 12,
                   column_split = factor(col_split, levels = unique(col_split))
)

# Show the heatmap
Heatmap <- draw(Heatmap)

svg("HH8_HH9_SOX8_Coexpression_Heatmap_RNA_assay_filtered.svg", width = 8, height = 6)
draw(Heatmap)
dev.off()

# add cluster labels and re-order
Correlation_matrix_filtered_and_annotated <- Extract_HM_Clusters(Heatmap, correlation_matrix_filtered_mat)

# Saving co-expression values together along with the cluster number and cell type annotations across stages
write.csv(Correlation_matrix_filtered_and_annotated, "Heatmap_clusters/Sox8_coexpression_clusters_HH8_HH9_JASPAR_RNA_assay_mean+stdv_12C_SOX8_cells.csv")
