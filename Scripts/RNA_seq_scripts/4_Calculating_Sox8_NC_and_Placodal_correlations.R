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

# Load and check the Seurat staged objects for NC_and_Placodal cells. Using those with Transcription factors only
HH7_NC_and_Placodal_TFs <- readRDS("HH7_NC_and_Placodal_TFs")
HH8_NC_and_Placodal_TFs <- readRDS("HH8_NC_and_Placodal_TFs")
HH9_NC_and_Placodal_TFs <- readRDS("HH9_NC_and_Placodal_TFs")
HH12_NC_and_Placodal_TFs <- readRDS("HH12_NC_and_Placodal_TFs")
HH14_NC_and_Placodal_TFs <- readRDS("HH14_NC_and_Placodal_TFs")
HH16_NC_and_Placodal_TFs <- readRDS("HH16_NC_and_Placodal_TFs")

# extract count data from NC_and_Placodal_only SEURAT object, and then subset based on TF genes only
HH7_NC_and_Placodal_TFs_counts <- t(HH7_NC_and_Placodal_TFs[["SCT"]]$counts)
HH8_NC_and_Placodal_TFs_counts <- t(HH8_NC_and_Placodal_TFs[["SCT"]]$counts)
HH9_NC_and_Placodal_TFs_counts <- t(HH9_NC_and_Placodal_TFs[["SCT"]]$counts)
HH12_NC_and_Placodal_TFs_counts <- t(HH12_NC_and_Placodal_TFs[["SCT"]]$counts)
HH14_NC_and_Placodal_TFs_counts <- t(HH14_NC_and_Placodal_TFs[["SCT"]]$counts)
HH16_NC_and_Placodal_TFs_counts <- t(HH16_NC_and_Placodal_TFs[["SCT"]]$counts)

# calculate correlation matrix and sort from high to low
HH7_NC_and_Placodal_TF_correlation_matrix <- WGCNA::cor(HH7_NC_and_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH7_NC_and_Placodal <- data.frame(as.list(HH7_NC_and_Placodal_TF_correlation_matrix["SOX8", ]))

HH8_NC_and_Placodal_TF_correlation_matrix <- WGCNA::cor(HH8_NC_and_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH8_NC_and_Placodal <- data.frame(as.list(HH8_NC_and_Placodal_TF_correlation_matrix["SOX8", ]))

HH9_NC_and_Placodal_TF_correlation_matrix <- WGCNA::cor(HH9_NC_and_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH9_NC_and_Placodal <- data.frame(as.list(HH9_NC_and_Placodal_TF_correlation_matrix["SOX8", ]))

HH12_NC_and_Placodal_TF_correlation_matrix <- WGCNA::cor(HH12_NC_and_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH12_NC_and_Placodal <- data.frame(as.list(HH12_NC_and_Placodal_TF_correlation_matrix["SOX8", ]))

HH14_NC_and_Placodal_TF_correlation_matrix <- WGCNA::cor(HH14_NC_and_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH14_NC_and_Placodal <- data.frame(as.list(HH14_NC_and_Placodal_TF_correlation_matrix["SOX8", ]))

HH16_NC_and_Placodal_TF_correlation_matrix <- WGCNA::cor(HH16_NC_and_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH16_NC_and_Placodal <- data.frame(as.list(HH16_NC_and_Placodal_TF_correlation_matrix["SOX8", ]))

# Make a dataframe from Sox8 Spearman's rank values for each stage
# First, assign rownames to each Sox8 correlation dataframe.
rownames(Sox8_co_expressed_TFs_HH7_NC_and_Placodal) <- "HH7"
rownames(Sox8_co_expressed_TFs_HH8_NC_and_Placodal) <- "HH8"
rownames(Sox8_co_expressed_TFs_HH9_NC_and_Placodal) <- "HH9"
rownames(Sox8_co_expressed_TFs_HH12_NC_and_Placodal) <- "HH12"
rownames(Sox8_co_expressed_TFs_HH14_NC_and_Placodal) <- "HH14"
rownames(Sox8_co_expressed_TFs_HH16_NC_and_Placodal) <- "HH16"
# Then combine the dataframes into 1 using rbind()
Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages <- rbind(Sox8_co_expressed_TFs_HH7_NC_and_Placodal,
                                                  Sox8_co_expressed_TFs_HH8_NC_and_Placodal,
                                                  Sox8_co_expressed_TFs_HH9_NC_and_Placodal,
                                                  Sox8_co_expressed_TFs_HH12_NC_and_Placodal,
                                                  Sox8_co_expressed_TFs_HH14_NC_and_Placodal,
                                                  Sox8_co_expressed_TFs_HH16_NC_and_Placodal)
Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages <- t(Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages)

# Remove SOX8 from dataframe and Replace all NA values with a 0
Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages <- subset(Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages, rownames(Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages) != "SOX8")
Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages[is.na(Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages)] <- 0

write.csv(Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages, "Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages.csv")

# Create the heatmap
pheatmap(Sox8_NC_and_Placodal_co_expressed_TFs_combined_stages,
         color = colorRampPalette(c("blue", "white", "red"))(25), # Blue to white to red
         scale = "none",  # Data is already scaled between -1 and 1
         cluster_rows = TRUE,  # Do not cluster rows, display as is
         cluster_cols = FALSE,  # Do not cluster columns, display as is
         display_numbers = FALSE,  # Optionally display the correlation values
         show_rownames = FALSE,
         main = "NC_and_Placodal Sox8 gene Correlations")  # Title of the heatmap


Sox8_co_expressed_TFs_HH7_NC_and_Placodal_sorted <- data.frame(as.list(sort(Sox8_co_expressed_TFs_HH7_NC_and_Placodal, decreasing = TRUE)))
