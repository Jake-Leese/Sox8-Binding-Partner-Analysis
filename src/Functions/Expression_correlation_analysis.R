# Large overall function for extracting correlation values, thresholding and visualising counts against correlations

# Load sub-functions 
source("/data/Sox8_binding_partner_analysis/src/Functions/subset_seurat_features.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_count_data.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Correlation_analysis.R")

# Function arguments:
# - Seurat_object = the starting seurat object
# - Gene_subset (optional) = a set of genes from which to use for the correlation analysis
#                            Can be given as a character vector or one of two specified TF databases, "JASPAR2024" or "AnimalTFDB"
#                            If not specified, full seurat feature list will be used
# - Stage = Stage by which to subset the cells used in the correaltion analysis
#           e.g. "HH7"
# - Cell_Type = Cell type by which to subset the cells used in the correlation analysis, based on "ectoderm_type" metadata
#               e.g. "placode"
# - Assay = count assay e.g. "RNA" or "SCT"
# - correlation_test = "spearman" or "pearson" 
# - Primary_Gene = specify the gene that you want to obtain correlation values against, e.g. "SOX8"
# - Subset_Cells_by_Gene = logical argument whether to use all the cells specified, or just those that have at least 1 count of the primary gene
# - threshold (optional) = Set a threshold for what should be considered co-expressed. 
#                          Can be given as a number < 0.0 - 1.0 > or "mean + stdev"
#                          If left, no threshold will be set
# - genes_to_label = a character vector containing genes of interest that you wish to label in the counts against correlation plots


Expression_correlation_analysis <- function(Seurat_obj, Gene_subset, Stage, Cell_Type, Assay, correlation_test, Primary_Gene, Subset_Cells_by_Gene, threshold, genes_to_label) {
  
  #
  if (missing(Gene_subset)) {
    Seurat <- Seurat_obj
  } else if (Gene_subset == "JASPAR2024") {
    Seurat <- subset_seurat_features(Seurat_obj, DB = "JASPAR2024")
  } else if (Gene_subset == "AnimalTFDB") {
    Seurat <- subset_seurat_features(Seurat_obj, DB = "AnimalTFDB")
  } else {
    Seurat <- DietSeurat(Seurat_obj, features = name_vector)
  }
  
  # Extract count data for specified stage and cell_type using custom function
  if (Subset_Cells_by_Gene == TRUE) {
    Counts <- Extract_count_data(Seurat, Stage = Stage, Cell_type = Cell_Type, Assay = Assay, data_type = "counts", subset_cells_by_gene = Primary_Gene)
  } else {
    Counts <- Extract_count_data(Seurat, Stage = Stage, Cell_type = Cell_Type, Assay = Assay, data_type = "counts")
  }
  
  # Extract correlation values against the specified primary gene using a custom function
  Gene_correlations <- Correlation_analysis(Counts, remove_null_count_genes = TRUE, Test = correlation_test, Gene = Primary_Gene)
  
  #################################################################
  #### Extract the counts and correlation values for each gene ####
  
  # Extract gene names from correlation_values
  Genes <- colnames(Gene_correlations)
  
  # Initialize an empty data frame to store the results
  result_df <- data.frame(matrix(ncol = length(Genes), nrow = 1))
  colnames(result_df) <- Genes
  
  
  for (name in Genes) {
    # Check if the column exists in counts
    if (name %in% colnames(Counts)) {
      # Sum the values under the corresponding column in counts
      column_sum <- sum(Counts[, name])
      
      # Add the sum to the new dataframe under the appropriate column
      result_df[1, name] <- column_sum
    } else {
      # If the column name is not found, you can set it to NA
      result_df[1, name] <- NA
    }
  }
  
  result_df <- t(result_df)
  correlation_vals <- as.numeric(Gene_correlations[1,])
  
  result_df <- as.data.frame(cbind(result_df, correlations = correlation_vals))
  colnames(result_df) <- c("total_counts", "correlations")
  result_df$log_total_counts <- log(result_df$total_counts)
  
  #####################################################################
  #### Add threshold and label columns based on function arguments ####
  
  if (missing(threshold)) {
    result_df <- result_df
  } else if (threshold == "mean + stdev") {
    result_df$threshold <- ifelse(result_df$correlations > (sd(result_df$correlations[result_df$correlations != 1]) + mean(result_df$correlations[result_df$correlations != 1])), "TRUE", "FALSE")
    message(sprintf("mean + stdev = %.5f", sd(result_df$correlations[result_df$correlations != 1]) + mean(result_df$correlations[result_df$correlations != 1])))
  } else {
    result_df$threshold <- ifelse(result_df$correlations > as.numeric(threshold), "TRUE", "FALSE")
  }
  
  if (missing(genes_to_label)) {
    result_df <- result_df
  } else {
    result_df$label <- ifelse(rownames(result_df) %in% genes_to_label, rownames(result_df), NA) # Adds a label column expressing TRUE/FALSE for the specified genes in the character vector "genes_to_label"
  }
    
  return(result_df)
  
}
