# Script to calculate correlation coefficients and write to files

setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(TFBSTools)
library(JASPAR2024)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(universalmotif)
library(ggplot2)

source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_HM_Clusters.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/subset_seurat_features.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_count_data.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Correlation_analysis.R")

# Load the original Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

# Subset Seurat object using custom subset_seurat_features() function that takes arguments
# Seurat_obj, DB ("JASPAR2024" (vertebrates) or "AnimalTFDB" (Ggal)) 

HHall_ectoderm_TFs <- subset_seurat_features(HHall_ectoderm, DB = "JASPAR2024")

#### EXTRACT COUNT DATA ####

# Custom function to extract count data for: Stage, Cell_type (ectoderm type), Sample (origin.ident)(optional), and Assay ("RNA" or "SCT"), data_type ("counts" or "data" (log.normalised)), and subset_cells_by_gene ("gene.name") 
NC_HH9_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "NC", Assay = "RNA", data_type = "counts")
NC_HH9_SCT_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "NC", Assay = "SCT", data_type = "counts")

dim(NC_HH9_SCT_Counts)
dim(NC_HH9_RNA_Counts)


#### CALCULATE CORRELATION COEFFICIENTS ####

# Custom function to calculate correlation coefficients 
# Function takes arguments: remove_null_count_genes (logical TRUE or FALSE), Test ("spearman" or "pearson"), Gene (Specify a gene to take all values against. If left blank, full correlation matrix will be returned)

# Returns whole correlation matrix
NC_HH9_RNA_pearson <- Correlation_analysis(NC_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "pearson")
NC_HH9_SCT_pearson <- Correlation_analysis(NC_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "pearson")
NC_HH9_RNA_spearman <- Correlation_analysis(NC_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman")
NC_HH9_SCT_spearman <- Correlation_analysis(NC_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "spearman")

# Returns correlation values for specified gene only
NC_HH9_RNA_pearson_SOX8 <- Correlation_analysis(NC_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "pearson", Gene = "SOX8")
NC_HH9_SCT_pearson_SOX8 <- Correlation_analysis(NC_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "pearson", Gene = "SOX8")
NC_HH9_RNA_spearman_SOX8 <- Correlation_analysis(NC_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
NC_HH9_SCT_spearman_SOX8 <- Correlation_analysis(NC_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")

# Retrieve stats about the correlation matrix
dim(NC_HH9_SCT_spearman_SOX8) # number of genes
median(abs(NC_HH9_RNA_pearson_SOX8)) # median abs value
max(NC_HH9_RNA_pearson_SOX8[NC_HH9_RNA_pearson_SOX8 < 0.99]) # Maximum correlation value
sd(NC_HH9_RNA_pearson_SOX8[NC_HH9_RNA_pearson_SOX8 < 0.99]) # standard deviation

# Retrieve stats about the correlation dataframe
median(abs(as.numeric(NC_HH9_SCT_pearson_SOX8[1,]))) # median abs value
max(as.numeric(NC_HH9_SCT_pearson_SOX8[NC_HH9_SCT_pearson_SOX8 < 0.99])) # Maximum correlation value
sd(as.numeric(NC_HH9_SCT_pearson_SOX8[NC_HH9_SCT_pearson_SOX8 < 0.99])) # standard deviation

#### PLOT CORRELATION COEFFICIENT DISTRIBUTIONS #### 

# extract values excl diagonals and any that are > 0.99 and return to a dataframe with correlation and method columns
NC_HH9_RNA_pearson_cor_val <- data.frame(correlation = (NC_HH9_RNA_pearson[upper.tri(NC_HH9_RNA_pearson)])[NC_HH9_RNA_pearson[upper.tri(NC_HH9_RNA_pearson)] < 0.99], method = "RNA_pearson")
NC_HH9_SCT_pearson_cor_val <- data.frame(correlation = (NC_HH9_SCT_pearson[upper.tri(NC_HH9_SCT_pearson)])[NC_HH9_SCT_pearson[upper.tri(NC_HH9_SCT_pearson)] < 0.99], method = "SCT_pearson")
NC_HH9_RNA_spearman_cor_val <- data.frame(correlation = (NC_HH9_RNA_spearman[upper.tri(NC_HH9_RNA_spearman)])[NC_HH9_RNA_spearman[upper.tri(NC_HH9_RNA_spearman)] < 0.99], method = "RNA_spearman")
NC_HH9_SCT_spearman_cor_val <- data.frame(correlation = (NC_HH9_SCT_spearman[upper.tri(NC_HH9_SCT_spearman)])[NC_HH9_SCT_spearman[upper.tri(NC_HH9_SCT_spearman)] < 0.99], method = "SCT_spearman")

# extract values from gene specific correlation dataframe and format appropriately for plotting
NC_HH9_RNA_pearson_cor_val <- data.frame(correlation = as.numeric(NC_HH9_RNA_pearson_SOX8[NC_HH9_RNA_pearson_SOX8 < 0.99]), method = "RNA_pearson") 
NC_HH9_SCT_pearson_cor_val <- data.frame(correlation = as.numeric(NC_HH9_SCT_pearson_SOX8[NC_HH9_SCT_pearson_SOX8 < 0.99]), method = "SCT_pearson")
NC_HH9_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(NC_HH9_RNA_spearman_SOX8[NC_HH9_RNA_spearman_SOX8 < 0.99]), method = "RNA_spearman")
NC_HH9_SCT_spearman_cor_val <- data.frame(correlation = as.numeric(NC_HH9_SCT_spearman_SOX8[NC_HH9_SCT_spearman_SOX8 < 0.99]), method = "SCT_spearman")

# combine dataframes 
combined_df <- rbind(NC_HH9_RNA_spearman_cor_val, NC_HH9_SCT_spearman_cor_val, NC_HH9_RNA_pearson_cor_val, NC_HH9_SCT_pearson_cor_val)

# Plot as box plots for each correlation matrix
ggplot(combined_df, aes(x = method, y = correlation)) + 
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) + # Box plot appearance
  geom_jitter(width = 0.3, color = "blue", size = 0.1, alpha = 0.7) +  # Add points for individual values
  labs(title = "HH9 NC cells SOX8 correlation coefficients", 
       x = "Method", y = "Correlation") + 
  theme_minimal()


write.csv(combined_df, "/data/Sox8_binding_partner_analysis/Plots/QC/SOX8+_cells_correlation_coefficient_vals.csv")
