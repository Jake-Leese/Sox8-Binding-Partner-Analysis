# Script to calculate correlation coefficients, plot and write them to file
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
placode_HH8_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH8", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
NC_HH8_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH8", Cell_type = "NC", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
placode_HH9_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
NC_HH9_RNA_Counts <- Extract_count_data(HHall_ectoderm_TFs, Stage = "HH9", Cell_type = "NC", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")

#dim(placode_HH9_SCT_Counts)
dim(placode_HH9_RNA_Counts)


#### CALCULATE CORRELATION COEFFICIENTS ####

# Custom function to calculate correlation coefficients 
# Function takes arguments: remove_null_count_genes (logical TRUE or FALSE), Test ("spearman" or "pearson"), Gene (Specify a gene to take all values against. If left blank, full correlation matrix will be returned)

# Returns whole correlation matrix
#placode_HH9_RNA_pearson <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "pearson")
#placode_HH9_SCT_pearson <- Correlation_analysis(placode_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "pearson")
#placode_HH9_RNA_spearman <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman")
#placode_HH9_SCT_spearman <- Correlation_analysis(placode_HH9_SCT_Counts, remove_null_count_genes = TRUE, Test = "spearman")

# Returns correlation values for specified gene only
placode_HH8_RNA_spearman_SOX8 <- Correlation_analysis(placode_HH8_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
NC_HH8_RNA_spearman_SOX8 <- Correlation_analysis(NC_HH8_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
placode_HH9_RNA_spearman_SOX8 <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
NC_HH9_RNA_spearman_SOX8 <- Correlation_analysis(NC_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")

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
#placode_HH9_RNA_spearman_cor_val <- data.frame(correlation = (placode_HH9_RNA_spearman[upper.tri(placode_HH9_RNA_spearman)])[placode_HH9_RNA_spearman[upper.tri(placode_HH9_RNA_spearman)] < 0.99], method = "RNA_spearman")
#placode_HH9_SCT_spearman_cor_val <- data.frame(correlation = (placode_HH9_SCT_spearman[upper.tri(placode_HH9_SCT_spearman)])[placode_HH9_SCT_spearman[upper.tri(placode_HH9_SCT_spearman)] < 0.99], method = "SCT_spearman")

# extract values from gene specific correlation dataframe and format appropriately for plotting
placode_HH8_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(placode_HH8_RNA_spearman_SOX8[placode_HH8_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH8_placode") 
NC_HH8_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(NC_HH8_RNA_spearman_SOX8[NC_HH8_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH8_NC")
placode_HH9_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(placode_HH9_RNA_spearman_SOX8[placode_HH9_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH9_placode")
NC_HH9_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(NC_HH9_RNA_spearman_SOX8[NC_HH9_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH9_NC")

# combine dataframes 
Combined_df <- rbind(placode_HH8_RNA_spearman_cor_val, NC_HH8_RNA_spearman_cor_val, placode_HH9_RNA_spearman_cor_val, NC_HH9_RNA_spearman_cor_val)

# Plot as box plots for each correlation matrix
ggplot(Combined_df, aes(x = Stage_CellType, y = correlation)) + 
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) + # Box plot appearance
  geom_jitter(width = 0.3, color = "blue", size = 0.1, alpha = 0.7) +  # Add points for individual values
  labs(title = "SOX8 RNA spearman correlations for SOX8+ cells", 
       x = "Method", y = "Correlation") + 
  theme_minimal()

# Calculate the mean + SD for each group
thresholds <- Combined_df %>%
  group_by(Stage_CellType) %>%
  summarize(threshold = mean(correlation) + sd(correlation))

# Join the threshold back to the original dataframe for comparison
Combined_df <- Combined_df %>%
  left_join(thresholds, by = "Stage_CellType") %>%
  mutate(above_threshold = factor(correlation > threshold, levels = c(FALSE, TRUE)))

# Plot correlation values as box plots and color points above the mean + stdev threshold as red vs blue.
plot <- ggplot(Combined_df, aes(x = Stage_CellType, y = correlation)) + 
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) + # Box plot appearance
  geom_jitter(aes(color = above_threshold), 
              width = 0.3, size = 0.1, alpha = 0.7) +  # Add points for individual values
  scale_color_manual(values = c("FALSE" = "darkblue", "TRUE" = "red")) +
  labs(title = "SOX8 correlation coefficients - RNA spearman SOX8+ cells", 
       x = "Stage_CellType", y = "Correlation") + 
  theme_minimal()

svg("/data/Sox8_binding_partner_analysis/Plots/QC/SOX8_correlation_coefficients_RNA_Sox8_positive.svg", width = 8, height = 6)
print(plot)
dev.off()

write.csv(Combined_df, "/data/Sox8_binding_partner_analysis/Plots/QC/SOX8_correlation_coefficients_HH9_placode_Sox8_positive.csv")


#### PLOTTING CORRELATION AGAINST TOTAL COUNTS FOR EACH GENE IN RELATION TO SOX8 ####

# custom function to combine total_counts and correlations for each gene based on a provided counts matrix and correlation df, respectively. 
# The function returns a dataframe with the gene, counts, correlations, and log_counts, and also returns a scatter plot of log total counts against correlation

HH8_RNA_spearman_placode_counts_and_cor <- plot_correlation_against_counts(placode_HH8_RNA_Counts, placode_HH8_RNA_spearman_SOX8)
HH8_RNA_spearman_NC_counts_and_cor <- plot_correlation_against_counts(NC_HH8_RNA_Counts, NC_HH8_RNA_spearman_SOX8)
HH9_RNA_spearman_placode_counts_and_cor <- plot_correlation_against_counts(placode_HH9_RNA_Counts, placode_HH9_RNA_spearman_SOX8)
HH9_RNA_spearman_NC_counts_and_cor <- plot_correlation_against_counts(NC_HH9_RNA_Counts, NC_HH9_RNA_spearman_SOX8)

# Comparing total SOX8 counts
log(sum(placode_HH8_RNA_Counts[, "SOX8"]))  # 6.338594
log(sum(NC_HH8_RNA_Counts[, "SOX8"]))   # 6.948897
log(sum(placode_HH9_RNA_Counts[, "SOX8"]))  # 6.472346 
log(sum(NC_HH9_RNA_Counts[, "SOX8"]))  # 9.013717

#write.csv(HH9_RNA_spearman_placode_counts_and_cor, "/data/Sox8_binding_partner_analysis/Plots/QC/Counts_against_correlation_HH9_RNA_spearman_placode_SOX8_cells.csv")

HH9_RNA_spearman_placode_counts_and_cor %>% arrange(desc(correlations))

HH8_placode_points_to_label <- c("LMX1B", "SP8", "TCF12", "DLX6", "TFAP2C", "LMX1A", "NR6A1", "NOTO", "SP5", "BHLHE23", "ASCL1") # For placode plots
HH8_NC_points_to_label <- c("SOX9", "ETS1", "TFAP2B", "CREB3L1", "FOXD3", "SOX10", "TFAP2E", "TFAP2A", "NEUROG2", "MYC", "FLI1") # for NC plots
HH9_placode_points_to_label <- c("PAX2", "GBX2", "MYCN", "NR6A1", "SOX13", "SIX1", "TFAP2A", "LMX1A", "FOXG1", "SOX10", "RXRG") # For placode plots
HH9_NC_points_to_label <- c("TFAP2B", "ETS1", "ATF4", "MYC", "NR6A1", "SOX10", "TCF3", "SOX4", "TCF12", "TFDP1") # for NC plots

# Adding a label column as a logical vector stating TRUE for the specified genes above. This is then called into ggplot for labelling those points
HH8_RNA_spearman_placode_counts_and_cor$label <- ifelse(rownames(HH8_RNA_spearman_placode_counts_and_cor) %in% HH8_placode_points_to_label, rownames(HH8_RNA_spearman_placode_counts_and_cor), NA)
HH8_RNA_spearman_NC_counts_and_cor$label <- ifelse(rownames(HH8_RNA_spearman_NC_counts_and_cor) %in% HH8_NC_points_to_label, rownames(HH8_RNA_spearman_NC_counts_and_cor), NA)
HH9_RNA_spearman_placode_counts_and_cor$label <- ifelse(rownames(HH9_RNA_spearman_placode_counts_and_cor) %in% HH9_placode_points_to_label, rownames(HH9_RNA_spearman_placode_counts_and_cor), NA)
HH9_RNA_spearman_NC_counts_and_cor$label <- ifelse(rownames(HH9_RNA_spearman_NC_counts_and_cor) %in% HH9_NC_points_to_label, rownames(HH9_RNA_spearman_NC_counts_and_cor), NA)

# Adding a threshold label to dataframes (TRUE or FALSE) based on whether the values fall above the stdev+mean threshold
HH8_RNA_spearman_placode_counts_and_cor$threshold <- ifelse(HH8_RNA_spearman_placode_counts_and_cor$correlations > (sd(HH8_RNA_spearman_placode_counts_and_cor$correlations) + mean(HH8_RNA_spearman_placode_counts_and_cor$correlations)), "TRUE", "FALSE")
HH8_RNA_spearman_NC_counts_and_cor$threshold <- ifelse(HH8_RNA_spearman_NC_counts_and_cor$correlations > (sd(HH8_RNA_spearman_NC_counts_and_cor$correlations) + mean(HH8_RNA_spearman_NC_counts_and_cor$correlations)), "TRUE", "FALSE")
HH9_RNA_spearman_placode_counts_and_cor$threshold <- ifelse(HH9_RNA_spearman_placode_counts_and_cor$correlations > (sd(HH9_RNA_spearman_placode_counts_and_cor$correlations) + mean(HH9_RNA_spearman_placode_counts_and_cor$correlations)), "TRUE", "FALSE")
HH9_RNA_spearman_NC_counts_and_cor$threshold <- ifelse(HH9_RNA_spearman_NC_counts_and_cor$correlations > (sd(HH9_RNA_spearman_NC_counts_and_cor$correlations) + mean(HH9_RNA_spearman_NC_counts_and_cor$correlations)), "TRUE", "FALSE")


# Plot a scatter plot, label specified factors and colour points above the threshold in red
ggplot(HH9_RNA_spearman_placode_counts_and_cor, aes(x = log_total_counts, y = correlations)) +
  geom_point(aes(color = threshold)) +
  scale_color_manual(values = c("FALSE" = "darkblue", "TRUE" = "red")) +
  geom_label_repel(aes(label = label), color = "black") +
  labs(x = "log_total_counts", y = "correlation", title = "HH9_RNA_spearman_placode_counts_and_cor") +
  theme_minimal()

# Save dataframes
write.csv(HH9_RNA_spearman_NC_counts_and_cor, "/data/Sox8_binding_partner_analysis/Plots/QC/Counts_against_correlation_HH9_RNA_spearman_NC_SOX8_cells.csv")

# Order dataframes by column 
#HH8_RNA_spearman_NC_counts_and_cor %>% arrange(log_total_counts) 
