# Script to define co-expression of transcription factors with SOX8
# It performs the following steps:
# - Calculates common metrics for transcription factors across the entire scRNA-seq dataset
# - Subsets seurat object by stage (HH8/HH9), cell type (NC/PPR), and expression of SOX8 (>1 SCT counts)
# - Calculates equivalent metrics for each subset to compare against the entire dataset
# - Calculates Log fold change for each transcription factor and filters those that are not co-expressed in at least one subset
# - Plots a heatmap of co-expressed genes and uses hierarchical clustering to group them based on their expression patterns
#
# Performed using personal ArchR_Seurat_R_4.4.1.sif container

setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(dplyr)
library(ggplot2)
library(JASPAR2024)
library(TFBSTools)
library(ggrepel)
library(MAST)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(WGCNA)
library(stringr)


# Load any custom functions
source("/data/Sox8_binding_partner_analysis/src/Functions/subset_seurat_features.R")

# Load Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")  

# Subset for transcription factors
HHall_ectoderm_TFs <- subset_seurat_features(HHall_ectoderm, DB = "JASPAR2024")
HH8_ectoderm_TFs <- subset(HHall_ectoderm_TFs, subset = stage == "HH8")
HH9_ectoderm_TFs <- subset(HHall_ectoderm_TFs, subset = stage == "HH9")

# Assign labels based on ectoderm type, stage, and SOX8 expression
HHall_ectoderm_TFs$cell_type_stage <- NA  # Create a new column

# Define cell identities based on ectoderm type, stage, and SOX8 expression
HHall_ectoderm_TFs$cell_type_stage[HHall_ectoderm_TFs$ectoderm_type == "placode" & 
                                     HHall_ectoderm_TFs$stage == "HH8" & 
                             GetAssayData(HHall_ectoderm_TFs, assay = "SCT", layer = "counts")["SOX8", ] > 0] <- "HH8_PPR_SOX8"

HHall_ectoderm_TFs$cell_type_stage[HHall_ectoderm_TFs$ectoderm_type == "NC" & 
                                     HHall_ectoderm_TFs$stage == "HH8" & 
                             GetAssayData(HHall_ectoderm_TFs, assay = "SCT", layer = "counts")["SOX8", ] > 0] <- "HH8_NC_SOX8"

HHall_ectoderm_TFs$cell_type_stage[HHall_ectoderm_TFs$ectoderm_type == "placode" & 
                                     HHall_ectoderm_TFs$stage == "HH9" & 
                             GetAssayData(HHall_ectoderm_TFs, assay = "SCT", layer = "counts")["SOX8", ] > 0] <- "HH9_PPR_SOX8"

HHall_ectoderm_TFs$cell_type_stage[HHall_ectoderm_TFs$ectoderm_type == "NC" & 
                                     HHall_ectoderm_TFs$stage == "HH9" & 
                             GetAssayData(HHall_ectoderm_TFs, assay = "SCT", layer = "counts")["SOX8", ] > 0] <- "HH9_NC_SOX8"

# Set new cell identities in Seurat object
Idents(HHall_ectoderm_TFs) <- HHall_ectoderm_TFs$cell_type_stage

# Check if labels are correctly assigned
table(Idents(HHall_ectoderm_TFs))

# Prepare SCT assay for marker detection
HHall_ectoderm_TFs <- PrepSCTFindMarkers(HHall_ectoderm_TFs)

# Calculate pairwise differential gene expression between Sox8 positive PPR and NC cells at each stage
DE_HH8_Sox8_cells <- FindMarkers(HHall_ectoderm_TFs, 
                                 ident.1 = "HH8_NC_SOX8", 
                                 ident.2 = "HH8_PPR_SOX8", 
                                 test.use = "MAST", 
                                 logfc.threshold = 0.25)

DE_HH9_Sox8_cells <- FindMarkers(HHall_ectoderm_TFs, 
                                 ident.1 = "HH9_NC_SOX8", 
                                 ident.2 = "HH9_PPR_SOX8", 
                                 test.use = "MAST", 
                                 logfc.threshold = 0.25)

HH8_PPR_DE_TFs <- DE_HH8_Sox8_cells[DE_HH8_Sox8_cells$avg_log2FC < 0 & DE_HH8_Sox8_cells$p_val_adj < 0.05, ]
HH8_NC_DE_TFs <- DE_HH8_Sox8_cells[DE_HH8_Sox8_cells$avg_log2FC > 0 & DE_HH8_Sox8_cells$p_val_adj < 0.05, ]
HH9_PPR_DE_TFs <- DE_HH9_Sox8_cells[DE_HH9_Sox8_cells$avg_log2FC < 0 & DE_HH9_Sox8_cells$p_val_adj < 0.05, ]
HH9_NC_DE_TFs <- DE_HH9_Sox8_cells[DE_HH9_Sox8_cells$avg_log2FC > 0 & DE_HH9_Sox8_cells$p_val_adj < 0.05, ]

# Saving DE_TFs to meme

galgal6_gtf <- import("/data/Sox8_binding_partner_analysis/genome_files/Gallus_gallus.GRCg6a.97.gtf")
galgal6_gtf_df <- as.data.frame(galgal6_gtf)

# Function to return character vector of gene names, gene_ids, or gene_list, and return gene name and corresponding gene_ids
Add_gene_IDs <- function(gtf, gene_list) {
  filtered_gtf <- gtf %>%
    filter(str_detect(tolower(gene_name), tolower(paste(gene_list, collapse = "|"))) | gene_id %in% gene_list | transcript_id %in% gene_list) %>%
    select(gene_id, gene_name) %>%
    distinct(gene_id, .keep_all = TRUE)
}  

HH8_PPR_DE_TFs_ids <- Add_gene_IDs(galgal6_gtf_df, rownames(HH8_PPR_DE_TFs))
HH8_NC_DE_TFs_ids <- Add_gene_IDs(galgal6_gtf_df, rownames(HH8_NC_DE_TFs))
HH9_PPR_DE_TFs_ids <- Add_gene_IDs(galgal6_gtf_df, rownames(HH9_PPR_DE_TFs))
HH9_NC_DE_TFs_ids <- Add_gene_IDs(galgal6_gtf_df, rownames(HH9_NC_DE_TFs))

# Saving as .csv file
write.csv(HH8_PPR_DE_TFs_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_PPR_DE_TFs.csv")
write.csv(HH8_NC_DE_TFs_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_DE_TFs.csv")
write.csv(HH9_PPR_DE_TFs_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_PPR_DE_TFs.csv")
write.csv(HH9_NC_DE_TFs_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_NC_DE_TFs.csv")

# Saving as gene_lists using gene_ids for enhancer annotation
writeLines(HH8_PPR_DE_TFs_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_PPR_DE_TFs.txt")
writeLines(HH8_NC_DE_TFs_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_DE_TFs.txt")
writeLines(HH9_PPR_DE_TFs_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_PPR_DE_TFs.txt")
writeLines(HH9_NC_DE_TFs_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_NC_DE_TFs.txt")


# Characterise transcription factors that have >0.25 Logfc against all other cells


# Merge DE results from HH8 and HH9
all_DE_TFs <- unique(c(rownames(DE_HH8_Sox8_cells), rownames(DE_HH9_Sox8_cells)))

# Scale genes and ensure that all features are kept. Many of "all_DE_TFs" are missing from original seurat scale.data.
HHall_ectoderm_TFs_scaled <- ScaleData(HHall_ectoderm_TFs, assay = "SCT", verbose = FALSE, features = rownames(HHall_ectoderm_TFs))

# Extract scaled expression for selected TFs
tf_expression <- GetAssayData(HHall_ectoderm_TFs_scaled, assay = "SCT", slot = "scale.data")[all_DE_TFs, ]

# Reorder columns based on cell type groups
tf_expression <- tf_expression[, order(Idents(HHall_ectoderm_TFs_scaled))]

# Define color palette
heatmap_colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Define annotation colors for cell types
cell_type_colors <- c("HH8_PPR_SOX8" = "#E69F00", 
                      "HH8_NC_SOX8" = "#56B4E9", 
                      "HH9_PPR_SOX8" = "#009E73", 
                      "HH9_NC_SOX8" = "#F0E442")

# Create column annotation (cell type colors)
col_ha <- HeatmapAnnotation(
  CellType = Idents(HHall_ectoderm_TFs_scaled),
  col = list(CellType = cell_type_colors),
  annotation_legend_param = list(title = "Cell Type")
)

# Create heatmap
Heatmap(tf_expression,
        name = "Expression",
        col = heatmap_colors,  # Apply the custom color scale
        show_row_names = TRUE,  # Show gene names
        show_column_names = FALSE, # Hide individual cell names
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        top_annotation = col_ha,  # Add cell type annotations
        column_split = Idents(HHall_ectoderm_TFs_scaled),  # Split by cell type
        row_names_gp = gpar(fontsize = 10),
        column_title = "Transcription Factor Expression in Sox8+ Cells"
)



HHall_ectoderm_TFs_scaled <- ScaleData(HHall_ectoderm_TFs, assay = "SCT", verbose = FALSE, features = rownames(HHall_ectoderm_TFs))


setdiff(all_DE_TFs, rownames(tf_expression))
length(missing_TFs)




# Function to calculate gene expression metrics
# Log2_mean is based on formula from this article. https://divingintogeneticsandgenomics.com/post/do-you-really-understand-log2fold-change-in-single-cell-rnaseq-data/#disqus_thread
# Using (exp() - 1) + 1e-9) accounts for genes that have low overall values
# Tested median and median absolute distribution, but even if we remove all 0 counts, 90% of TFs have a median of 1
calculate_metrics <- function(seurat_obj) {
  data <- as.matrix(GetAssayData(seurat_obj, assay = "SCT", layer = "counts"))
  metrics <- data.frame(
    Gene = rownames(data),
    Total_Counts = rowSums(data),
    Mean_Expression = rowMeans(data),
    Std_Dev = apply(data, 1, sd),
    log2_mean = log2(rowMeans(exp(data) - 1) + 1e-9),
    Proportion_Expressing_Cells = rowMeans(data > 0)
  )
  return(metrics)
}
# Calculate metrics for the full dataset
full_metrics <- calculate_metrics(HHall_ectoderm_TFs)
HH8_metrics <- calculate_metrics(HH8_ectoderm_TFs)
HH9_metrics <- calculate_metrics(HH9_ectoderm_TFs)

# Subset based on SOX8 expression and cell type/stage
SOX8_positive <- HHall_ectoderm_TFs@assays$SCT@counts["SOX8", ] > 1  # Define SOX8+ cells

HH8_PPR_SOX8 <- subset(HHall_ectoderm_TFs, cells = WhichCells(HHall_ectoderm_TFs, expression = stage == "HH8" & ectoderm_type == "placode" & SOX8_positive))
HH8_NC_SOX8  <- subset(HHall_ectoderm_TFs, cells = WhichCells(HHall_ectoderm_TFs, expression = stage == "HH8" & ectoderm_type == "NC" & SOX8_positive))
HH9_PPR_SOX8 <- subset(HHall_ectoderm_TFs, cells = WhichCells(HHall_ectoderm_TFs, expression = stage == "HH9" & ectoderm_type == "placode" & SOX8_positive))
HH9_NC_SOX8  <- subset(HHall_ectoderm_TFs, cells = WhichCells(HHall_ectoderm_TFs, expression = stage == "HH9" & ectoderm_type == "NC" & SOX8_positive))

# Compute metrics for each subset
HH8_PPR_metrics <- calculate_metrics(HH8_PPR_SOX8)
HH8_NC_metrics  <- calculate_metrics(HH8_NC_SOX8)
HH9_PPR_metrics <- calculate_metrics(HH9_PPR_SOX8)
HH9_NC_metrics  <- calculate_metrics(HH9_NC_SOX8)

HH8_PPR_metrics %>% arrange(desc(Proportion_Expressing_Cells))
HH9_PPR_metrics %>% arrange(desc(Proportion_Expressing_Cells))

# Calculate log fold changes
calculate_logFC <- function(subset_metrics, full_metrics) {
  merged_data <- merge(subset_metrics, full_metrics, by = "Gene", suffixes = c("_subset", "_full"))
  merged_data$LogFC <- merged_data$log2_mean_subset - merged_data$log2_mean_full
  return(merged_data)
}

HH8_PPR_logFC <- calculate_logFC(HH8_PPR_metrics, HH8_metrics)
HH8_NC_logFC  <- calculate_logFC(HH8_NC_metrics, HH8_metrics)
HH9_PPR_logFC <- calculate_logFC(HH9_PPR_metrics, HH9_metrics)
HH9_NC_logFC  <- calculate_logFC(HH9_NC_metrics, HH9_metrics)

# Check top results for co-expression
HH8_PPR_logFC %>% arrange(desc(LogFC))
HH8_NC_logFC %>% arrange(desc(LogFC))
HH9_PPR_logFC %>% arrange(desc(LogFC))
HH9_NC_logFC %>% arrange(desc(LogFC)) 

# Merge stage/cell_type data
merged_LogFC_data <- data.frame(
  HH8_PPR = HH8_PPR_logFC$LogFC,
  HH8_NC  = HH8_NC_logFC$LogFC,
  HH9_PPR = HH9_PPR_logFC$LogFC,
  HH9_NC  = HH9_NC_logFC$LogFC
)

# Add gene names
rownames(merged_LogFC_data) <- HH8_PPR_logFC$Gene

# Filter: Keep genes with at least one LogFC > 0.5 in any of the four datasets
filtered_data <- merged_LogFC_data[apply(merged_LogFC_data, 1, function(x) any(x > 0.5, na.rm = TRUE)), ]

# Define co-expression threshold
threshold <- 0.5
Sox8_co_expressed <- unique(c(
  HH8_PPR_logFC$Gene[HH8_PPR_logFC$LogFC > threshold],
  HH8_NC_logFC$Gene[HH8_NC_logFC$LogFC > threshold],
  HH9_PPR_logFC$Gene[HH9_PPR_logFC$LogFC > threshold],
  HH9_NC_logFC$Gene[HH9_NC_logFC$LogFC > threshold]
))

# Convert to z-scores
expression_matrix <- as.matrix(GetAssayData(HHall_ectoderm_TFs, assay = "SCT", slot = "scale.data"))
z_scores <- t(apply(expression_matrix[Sox8_co_expressed, , drop = FALSE], 1, scale))

# Plot heatmap
pheatmap(z_scores, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = FALSE)
