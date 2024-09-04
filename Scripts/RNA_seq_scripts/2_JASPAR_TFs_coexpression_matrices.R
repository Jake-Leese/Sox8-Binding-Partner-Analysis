setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
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

head(HHall_ectoderm)

# Select the NC cells to be used for co-expression analysis
HHall_NC_only <- HHall_ectoderm[,HHall_ectoderm$ectoderm_type %in% "NC"]

head(HHall_NC_only)

# View processed and normalised count data for each cell
head(HHall_NC_only[["SCT"]]$counts)

### Subset count data for TFs only
# Download JASPAR2020 motif list for vertebrates
motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))
     
# Extract TF names from motif database
name_vector <- c()
for (i in 1:length(motifList)){
  name <- name(motifList[[i]])
  name_vector <- c(name_vector, name)
}

name_vector

# extract list of genes from Seurat object and check overlap with TF list
RNAseq_all_genes <- rownames(HHall_NC_only)
tf_genes <- intersect(RNAseq_all_genes, name_vector)

# extract count data from NC_only SEURAT object, and then subset based on TF genes only
NC_count_data <- HHall_NC_only[["SCT"]]$counts
NC_TF_count_data <- NC_count_data[rownames(NC_count_data) %in% tf_genes, ]
transformed_NC_TF_countr_data <- t(NC_TF_count_data)

# calculate correlation matrix
NC_TF_correlation_matrix <- WGCNA::cor(transformed_NC_TF_countr_data, method = "spearman")
Sox8_co_expressed_TFs_NC <- NC_TF_correlation_matrix["PCDH20", ]
Sox8_co_expression_Spearmans_sorted <- sort(Sox8_co_expressed_TFs_NC, decreasing = TRUE)
Sox8_co_expression_Spearmans_Positive <- data.frame(as.list(Sox8_co_expression_Spearmans_sorted[which(Sox8_co_expression_Spearmans_sorted > 0)]))
write.csv(Sox8_co_expression_Spearmans_Positive, "Sox8_co_expression_Spearmans_Positive.csv")



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
