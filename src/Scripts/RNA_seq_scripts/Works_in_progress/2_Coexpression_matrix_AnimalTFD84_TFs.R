# Base script to subset RNAseq dataset, filter TFs and create a correlation matrix based on their expression

setwd("/data/scRNAseq_objects")
.libPaths("/R/libs/R_scvi_integration_v3.4")

library(Seurat)
library(usethis)
library(devtools)
library(CSCORE)
library(TFBSTools)
library(JASPAR2020)
library(WGCNA)

# Load and check the Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

# Select the NC cells to be used for co-expression analysis
HHall_NC_only <- HHall_ectoderm[,HHall_ectoderm$ectoderm_type %in% "NC"]

# View processed and normalised count data for each cell
head(HHall_NC_only[["SCT"]]$counts)

# load AnimalTFDB TF database for Gallus gallus and extract gene (symbol) names as a vector
GalGal_TFs_df <- read.table("Gallus_gallus_TF.txt", sep = "\t", header = TRUE)
GalGal_TFs <- GalGal_TFs_df$Symbol

# extract list of genes from Seurat object and retrieve overlapping TFs. 814 out of the 1058 TFs from the animalTFD84 database. The ones
# that aren't mapped are most likely just those that are missing a name.
RNAseq_all_genes <- rownames(HHall_NC_only)
NC_only_TF_genes <- intersect(RNAseq_all_genes, GalGal_TFs)

# extract count data from NC_only SEURAT object, and then subset based on TF genes only
NC_count_data <- HHall_NC_only[["SCT"]]$counts
NC_TF_count_data <- NC_count_data[rownames(NC_count_data) %in% NC_only_TF_genes, ]
transformed_NC_TF_countr_data <- t(NC_TF_count_data)

# calculate correlation matrix
NC_TF_correlation_matrix <- WGCNA::cor(transformed_NC_TF_countr_data, method = "spearman")
Sox8_co_expressed_TFs_NC <- NC_TF_correlation_matrix["SOX8", ]
Sox8_co_expression_Spearmans_sorted <- sort(Sox8_co_expressed_TFs_NC, decreasing = TRUE)
Sox8_co_expression_Spearmans_Positive <- data.frame(as.list(Sox8_co_expression_Spearmans_sorted[which(Sox8_co_expression_Spearmans_sorted > 0)]))
write.csv(Sox8_co_expression_Spearmans_Positive, "Sox8_AnimalTFD84_co_expression_Spearmans_Positive.csv")



## plot data to heatmap
melted_correlation_matrix <- melt(NC_TF_correlation_matrix)

ggplot(data = melted_correlation_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 1, hjust = 1),
        axis.text.y = element_text(size = 1)) +
  coord_fixed()