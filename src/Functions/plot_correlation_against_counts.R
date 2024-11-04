# Function to make a scatter plot of correlation vs counts using correlation values calulated against one gene.

plot_correlation_against_counts <- function(count_matrix, correlation_values) {
  
  # Extract TF names from correlation_values
  TFs <- colnames(correlation_values)
  
  # Initialize an empty data frame to store the results
  result_df <- data.frame(matrix(ncol = length(TFs), nrow = 1))
  colnames(result_df) <- TFs
  
  
  for (name in TFs) {
    # Check if the column exists in counts
    if (name %in% colnames(count_matrix)) {
      # Sum the values under the corresponding column in counts
      column_sum <- sum(count_matrix[, name])
      
      # Add the sum to the new dataframe under the appropriate column
      result_df[1, name] <- column_sum
    } else {
      # If the column name is not found, you can set it to NA
      result_df[1, name] <- NA
    }
  }
  
  result_df <- t(result_df)
  correlation_vals <- as.numeric(correlation_values[1,])
  
  result_df <- as.data.frame(cbind(result_df, correlations = correlation_vals))
  colnames(result_df) <- c("total_counts", "correlations") 
  result_df <- subset(result_df, correlations != 1)
  result_df$log_total_counts <- log(result_df$total_counts)
  
  # Capture the name of the correlation_values dataframe as a string
  correlation_values_name <- deparse(substitute(correlation_values))
  
  # ggplot2 scatter plot
  plot <- ggplot(result_df, aes(x = log_total_counts, y = correlations)) +
    geom_point(color = "blue") +
    labs(x = "log_total_counts", y = "correlation", title = paste("Scatter Plot for", correlation_values_name)) +
    theme_minimal()
  
  print(plot)
  return(result_df)
}
