# Script to calculate correlation coefficients and write to files

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
placode_HH9_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
#placode_HH9_SCT_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "placode", Assay = "SCT", data_type = "counts", subset_cells_by_gene = "SOX8")

#dim(placode_HH9_SCT_Counts)
dim(placode_HH9_RNA_Counts)


#### CALCULATE CORRELATION COEFFICIENTS ####

# Custom function to calculate correlation coefficients 
# Function takes arguments: remove_null_count_genes (logical TRUE or FALSE), Test ("spearman" or "pearson"), Gene (Specify a gene to take all values against. If left blank, full correlation matrix will be returned)

# Returns whole correlation matrix
#placode_HH9_RNA_pearson <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "pearson")
#placode_HH9_SCT_pearson <- Correlation_analysis(placode_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "pearson")
placode_HH9_RNA_spearman <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman")
#placode_HH9_SCT_spearman <- Correlation_analysis(placode_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "spearman")

# Returns correlation values for specified gene only
#placode_HH9_RNA_pearson_SOX8 <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "pearson", Gene = "SOX8")
#placode_HH9_SCT_pearson_SOX8 <- Correlation_analysis(placode_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "pearson", Gene = "SOX8")
placode_HH9_RNA_spearman_SOX8 <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
#placode_HH9_SCT_spearman_SOX8 <- Correlation_analysis(placode_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")

# Retrieve stats about the correlation matrix
dim(placode_HH9_SCT_spearman_SOX8) # number of genes
median(abs(placode_HH9_RNA_pearson_SOX8)) # median abs value
max(placode_HH9_RNA_pearson_SOX8[placode_HH9_RNA_pearson_SOX8 < 0.99]) # Maximum correlation value
sd(placode_HH9_RNA_pearson_SOX8[placode_HH9_RNA_pearson_SOX8 < 0.99]) # standard deviation

# Retrieve stats about the correlation dataframe (specified gene)
dim(placode_HH9_RNA_pearson_SOX8)
median(abs(as.numeric(placode_HH9_SCT_pearson_SOX8[1,]))) # median abs value
max(as.numeric(placode_HH9_SCT_pearson_SOX8[placode_HH9_SCT_pearson_SOX8 < 0.99])) # Maximum correlation value
sd(as.numeric(placode_HH9_SCT_pearson_SOX8[placode_HH9_SCT_pearson_SOX8 < 0.99])) # standard deviation

#### PLOT CORRELATION COEFFICIENT DISTRIBUTIONS #### 

# extract values excl diagonals and any that are > 0.99 and return to a dataframe with correlation and method columns
#placode_HH9_RNA_pearson_cor_val <- data.frame(correlation = (placode_HH9_RNA_pearson[upper.tri(placode_HH9_RNA_pearson)])[placode_HH9_RNA_pearson[upper.tri(placode_HH9_RNA_pearson)] < 0.99], method = "RNA_pearson")
#placode_HH9_SCT_pearson_cor_val <- data.frame(correlation = (placode_HH9_SCT_pearson[upper.tri(placode_HH9_SCT_pearson)])[placode_HH9_SCT_pearson[upper.tri(placode_HH9_SCT_pearson)] < 0.99], method = "SCT_pearson")
placode_HH9_RNA_spearman_cor_val <- data.frame(correlation = (placode_HH9_RNA_spearman[upper.tri(placode_HH9_RNA_spearman)])[placode_HH9_RNA_spearman[upper.tri(placode_HH9_RNA_spearman)] < 0.99], method = "RNA_spearman")
#placode_HH9_SCT_spearman_cor_val <- data.frame(correlation = (placode_HH9_SCT_spearman[upper.tri(placode_HH9_SCT_spearman)])[placode_HH9_SCT_spearman[upper.tri(placode_HH9_SCT_spearman)] < 0.99], method = "SCT_spearman")

# extract values from gene specific correlation dataframe and format appropriately for plotting
#placode_HH9_RNA_pearson_cor_val <- data.frame(correlation = as.numeric(placode_HH9_RNA_pearson_SOX8[placode_HH9_RNA_pearson_SOX8 < 0.99]), method = "RNA_pearson") 
#placode_HH9_SCT_pearson_cor_val <- data.frame(correlation = as.numeric(placode_HH9_SCT_pearson_SOX8[placode_HH9_SCT_pearson_SOX8 < 0.99]), method = "SCT_pearson")
placode_HH9_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(placode_HH9_RNA_spearman_SOX8[placode_HH9_RNA_spearman_SOX8 < 0.99]), method = "RNA_spearman")
#placode_HH9_SCT_spearman_cor_val <- data.frame(correlation = as.numeric(placode_HH9_SCT_spearman_SOX8[placode_HH9_SCT_spearman_SOX8 < 0.99]), method = "SCT_spearman")

# combine dataframes 
combined_df <- rbind(placode_HH9_RNA_spearman_cor_val, placode_HH9_SCT_spearman_cor_val, placode_HH9_RNA_pearson_cor_val, placode_HH9_SCT_pearson_cor_val)

# Plot as box plots for each correlation matrix
ggplot(combined_df, aes(x = method, y = correlation)) + 
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) + # Box plot appearance
  geom_jitter(width = 0.3, color = "blue", size = 0.1, alpha = 0.7) +  # Add points for individual values
  labs(title = "HH9 placode cells SOX8 correlation coefficients", 
       x = "Method", y = "Correlation") + 
  theme_minimal()


write.csv(combined_df, "/data/Sox8_binding_partner_analysis/Plots/QC/SOX8_correlation_coefficients_HH9_placode_Sox8_positive.csv")


#### PLOTTING CORRELATION AGAINST TOTAL COUNTS FOR EACH GENE IN RELATION TO SOX8 ####

# custom function to combine total_counts and correlations for each gene based on a provided counts matrix and correlation df, respectively. 
# The function returns a dataframe with the gene, counts, correlations, and log_counts, and also returns a scatter plot of log total counts against correlation

HH9_RNA_spearman_placode_counts_and_cor <- plot_correlation_against_counts(placode_HH9_RNA_Counts, placode_HH9_RNA_spearman_SOX8)

log(sum(placode_HH9_RNA_Counts[, "SOX8"]))

write.csv(HH9_RNA_spearman_placode_counts_and_cor, "/data/Sox8_binding_partner_analysis/Plots/QC/Counts_against_correlation_HH9_RNA_spearman_placode_SOX8_cells.csv")

HH9_RNA_spearman_placode_counts_and_cor %>% arrange(desc(correlations))

# points_to_label <- c("PAX2", "GBX2", "MYCN", "NR6A1", "SOX13", "HES5", "ATF4", "LMX1A", "TCF12", "SOX21") # For placode plots
points_to_label <- c("TFAP2B", "ETS1", "ATF4", "MYC", "NR6A1", "SOX10", "TCF3", "SOX4", "TCF12", "TFDP1") # for NC plots

HH9_RNA_spearman_placode_counts_and_cor$label <- ifelse(rownames(HH9_RNA_spearman_placode_counts_and_cor) %in% points_to_label, rownames(HH9_RNA_spearman_placode_counts_and_cor), NA)

ggplot(HH9_RNA_spearman_placode_counts_and_cor, aes(x = log_total_counts, y = correlations)) +
  geom_point(color = "blue") +
  geom_label_repel(aes(label = label), color = "black") +
  labs(x = "log_total_counts", y = "correlation", title = "HH9_RNA_spearman_placode_counts_and_cor") +
  theme_minimal()

