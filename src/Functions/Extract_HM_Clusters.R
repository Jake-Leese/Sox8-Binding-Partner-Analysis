# Function to extract clusters from a heatmap 

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
  return(out)
}
