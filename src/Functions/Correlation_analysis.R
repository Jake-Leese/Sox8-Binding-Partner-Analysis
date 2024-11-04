# Function to calculate correlation coefficients from a scRNA-seq count matrix

Correlation_analysis <- function(count_matrix, remove_null_count_genes, Test, Gene) {
  
  if (remove_null_count_genes == TRUE) {
    count_matrix_processed <- count_matrix[,colSums(count_matrix) > 0]
  } else {
    count_matrix_processed <- count_matrix
  }
  
  correlation_matrix <- WGCNA::cor(count_matrix_processed, method = Test)
  
  if (missing(Gene)) {
    return(correlation_matrix)
  } else {
    Gene_correlation <- data.frame(as.list(correlation_matrix[Gene, ]))
    return(Gene_correlation)
  }
}
