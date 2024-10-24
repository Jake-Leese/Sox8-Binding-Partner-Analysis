setwd("/data/scRNAseq_objects")
.libPaths("/R/libs/R_scvi_integration_v3.4")

library(Seurat)
library(usethis)
library(devtools)
library(CSCORE)
library(TFBSTools)
library(JASPAR2020)
library(WGCNA)
library(pheatmap)

# Load and check the Seurat staged objects for NC cells. Using those with Transcription factors only
HH7_NC_TFs <- readRDS("HH7_NC_TFs")
HH8_NC_TFs <- readRDS("HH8_NC_TFs")
HH9_NC_TFs <- readRDS("HH9_NC_TFs")
HH12_NC_TFs <- readRDS("HH12_NC_TFs")
HH14_NC_TFs <- readRDS("HH14_NC_TFs")
HH16_NC_TFs <- readRDS("HH16_NC_TFs")

# extract count data from NC_only SEURAT object, and then subset based on TF genes only
HH7_NC_TFs_counts <- t(HH7_NC_TFs[["SCT"]]$counts)
HH8_NC_TFs_counts <- t(HH8_NC_TFs[["SCT"]]$counts)
HH9_NC_TFs_counts <- t(HH9_NC_TFs[["SCT"]]$counts)
HH12_NC_TFs_counts <- t(HH12_NC_TFs[["SCT"]]$counts)
HH14_NC_TFs_counts <- t(HH14_NC_TFs[["SCT"]]$counts)
HH16_NC_TFs_counts <- t(HH16_NC_TFs[["SCT"]]$counts)

# calculate correlation matrix and sort from high to low
HH7_NC_TF_correlation_matrix <- WGCNA::cor(HH7_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH7_NC <- data.frame(as.list(HH7_NC_TF_correlation_matrix["SOX8", ]))

HH8_NC_TF_correlation_matrix <- WGCNA::cor(HH8_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH8_NC <- data.frame(as.list(HH8_NC_TF_correlation_matrix["SOX8", ]))

HH9_NC_TF_correlation_matrix <- WGCNA::cor(HH9_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH9_NC <- data.frame(as.list(HH9_NC_TF_correlation_matrix["SOX8", ]))

HH12_NC_TF_correlation_matrix <- WGCNA::cor(HH12_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH12_NC <- data.frame(as.list(HH12_NC_TF_correlation_matrix["SOX8", ]))

HH14_NC_TF_correlation_matrix <- WGCNA::cor(HH14_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH14_NC <- data.frame(as.list(HH14_NC_TF_correlation_matrix["SOX8", ]))

HH16_NC_TF_correlation_matrix <- WGCNA::cor(HH16_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH16_NC <- data.frame(as.list(HH16_NC_TF_correlation_matrix["SOX8", ]))

# Make a dataframe from Sox8 Spearman's rank values for each stage
# First, assign rownames to each Sox8 correlation dataframe.
rownames(Sox8_co_expressed_TFs_HH7_NC) <- "HH7"
rownames(Sox8_co_expressed_TFs_HH8_NC) <- "HH8"
rownames(Sox8_co_expressed_TFs_HH9_NC) <- "HH9"
rownames(Sox8_co_expressed_TFs_HH12_NC) <- "HH12"
rownames(Sox8_co_expressed_TFs_HH14_NC) <- "HH14"
rownames(Sox8_co_expressed_TFs_HH16_NC) <- "HH16"
# Then combine the dataframes into 1 using rbind()
Sox8_NC_co_expressed_TFs_combined_stages <- rbind(Sox8_co_expressed_TFs_HH7_NC,
                                                  Sox8_co_expressed_TFs_HH8_NC,
                                                  Sox8_co_expressed_TFs_HH9_NC,
                                                  Sox8_co_expressed_TFs_HH12_NC,
                                                  Sox8_co_expressed_TFs_HH14_NC,
                                                  Sox8_co_expressed_TFs_HH16_NC)
Sox8_NC_co_expressed_TFs_combined_stages <- t(Sox8_NC_co_expressed_TFs_combined_stages)

# Remove SOX8 from dataframe and Replace all NA values with a 0
Sox8_NC_co_expressed_TFs_combined_stages <- subset(Sox8_NC_co_expressed_TFs_combined_stages, rownames(Sox8_NC_co_expressed_TFs_combined_stages) != "SOX8")
Sox8_NC_co_expressed_TFs_combined_stages[is.na(Sox8_NC_co_expressed_TFs_combined_stages)] <- 0

write.csv(Sox8_NC_co_expressed_TFs_combined_stages, "Sox8_NC_co_expressed_TFs_combined_stages.csv")

# Create the heatmap
pheatmap(Sox8_NC_co_expressed_TFs_HH7_HH8,
         color = colorRampPalette(c("blue", "white", "red"))(25), # Blue to white to red
         scale = "none",  # Data is already scaled between -1 and 1
         cluster_rows = TRUE,  # Do not cluster rows, display as is
         cluster_cols = FALSE,  # Do not cluster columns, display as is
         display_numbers = FALSE,  # Optionally display the correlation values
         fontsize_row = 2,
         main = "NC Sox8 gene Correlations")  # Title of the heatmap


Sox8_co_expressed_TFs_HH7_NC_sorted <- data.frame(as.list(sort(Sox8_co_expressed_TFs_HH7_NC, decreasing = TRUE)))
