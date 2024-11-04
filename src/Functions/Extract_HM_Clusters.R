# Function to extract clusters from a heatmap and assign them  

Extract_HM_Clusters <- function(Heatmap, correlation_matrix) {
  
  # Extract clusters and gene numbers 
  Cluster_list <- row_order(Heatmap)
  
  # Loop to extract corresponding gene names for each cluster
  for (i in 1:length(Cluster_list)) {
    if (i == 1) {                                                   # This if statement provides one iteration for the initial out object to be created without overwriting the updated version each time
      clu <- t(t(row.names(correlation_matrix)[Cluster_list[[i]]])) # extracts gene names that correspond to the row numbers of cluster, i.
      out <- cbind(clu, paste("Cluster", i, sep="_"))
      colnames(out) <- c("Gene", "Cluster")
    } else {
      clu <- t(t(row.names(correlation_matrix)[Cluster_list[[i]]]))
      clu <- cbind(clu, paste("Cluster", i, sep="_"))
      out <- rbind(out, clu)
    }
  }
  
  # Add a cluster label column to the original correlation_matrix
  correlation_matrix_clust <- as.data.frame(out[match(rownames(correlation_matrix),out[,1]), ]) # 
  correlation_matrix_df <- as.data.frame(correlation_matrix)
  correlation_matrix_df$Heatmap_Cluster <- c(correlation_matrix_clust[,2])
  
  # re-ordering based on cluster number
  correlation_matrix_df$Heatmap_Cluster <- as.numeric(sub(".*_(\\d+)$", "\\1", correlation_matrix_df$Heatmap_Cluster)) # Changing the cluster column to contain only numbers. This allows us to re-order the dataframe based on the order of the clusters
  correlation_matrix_df <- correlation_matrix_df[order(correlation_matrix_df$Heatmap_Cluster),]
  
  return(correlation_matrix_df)
}
