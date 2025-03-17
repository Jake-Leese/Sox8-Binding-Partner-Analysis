# Simplified script to calculate correlation coefficients and visualise distributions using custom function, Expression_correlation_analysis
# Run using ArchR_Seurat_R_4.4.1.sif container

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

source("/data/Sox8_binding_partner_analysis/src/Functions/Expression_correlation_analysis.R")

# Load Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

# Subset to transcription factors only (optional)
HHall_ectoderm_TFs <- subset_seurat_features(HHall_ectoderm, DB = "JASPAR2024")

# Specify TFs of interest to label in plot
HH8_placode_points_to_label <- c("LMX1B", "SP8", "TCF12", "DLX6", "TFAP2C", "LMX1A", "NR6A1", "NOTO", "SP5", "BHLHE23", "ASCL1", "PRDM1")
HH8_NC_points_to_label <- c("SOX9", "ETS1", "TFAP2B", "CREB3L1", "FOXD3", "SOX10", "TFAP2E", "TFAP2A", "NEUROG2", "MYC", "FLI1") # for NC plots
HH9_placode_points_to_label <- c("PAX2", "GBX2", "MYCN", "NR6A1", "SOX13", "HES5", "ATF4", "LMX1A", "TCF12", "SOX21", "PRDM1") # For placode plots
HH9_NC_points_to_label <- c("TFAP2B", "ETS1", "ATF4", "MYC", "NR6A1", "SOX10", "TCF3", "SOX4", "TCF12", "TFDP1") # for NC plots

HH8_Placode_SOX8_counts_v_correlation_TFs <- Expression_correlation_analysis(Seurat_obj = HHall_ectoderm_TFs,
                                                                         Stage = "HH8",
                                                                         Cell_Type = "placode",
                                                                         Assay = "RNA",
                                                                         correlation_test = "spearman",
                                                                         Primary_Gene = "SOX8",
                                                                         Subset_Cells_by_Gene = FALSE,
                                                                         threshold = "mean + stdev",
                                                                         genes_to_label = HH8_placode_points_to_label)



# Remove the primary_gene value (=1) before plotting
HH8_Placode_SOX8_counts_v_correlation_TFs_plot <- HH8_Placode_SOX8_counts_v_correlation_TFs[HH8_Placode_SOX8_counts_v_correlation_TFs$correlations != 1, ]

# Plot a scatter plot, label specified factors and colour points above the threshold in red
ggplot(HH8_Placode_SOX8_counts_v_correlation_TFs_plot, aes(x = log_total_counts, y = correlations)) +
  geom_point(aes(color = threshold)) +
  scale_color_manual(values = c("FALSE" = "darkblue", "TRUE" = "red")) +
  geom_label_repel(aes(label = label), color = "black") +
  labs(x = "log_total_counts", y = "correlation", title = "HH8_Placode_SOX8_counts_v_correlation_TFs_plot") +
  theme_minimal()


# Extract list of genes above threshold
HH9_Placode_SOX8_co_expressed_all_genes <- HH9_Placode_SOX8_counts_v_correlation_all_genes[HH9_Placode_SOX8_counts_v_correlation_all_genes$threshold == TRUE, ] 
HH9_Placode_SOX8_co_expressed_all_genes <- rownames(HH9_Placode_SOX8_co_expressed_all_genes)

# Save character vector as a .txt file, where each string is given per line
writeLines(HH9_Placode_SOX8_co_expressed_all_genes, "/data/Sox8_binding_partner_analysis/Enhancer_annotation/Gene_lists/HH9_Placode_SOX8_co-expressed_mean_stdev_0.159_all_genes.txt")
write.tx
