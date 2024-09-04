# Function to extract peak matrix from an ArchR project and label rows by peak names

ExtractArchRPeakMatrix <- function(ArchR_project) {
  
  Full_peakMatrix <- getMatrixFromProject(ArchR_project, useMatrix = "PeakMatrix")      
  Full_peakMatrix_reads <- Full_peakMatrix@assays@data$PeakMatrix  
  Full_peakMatrix_names <- Full_peakMatrix@rowRanges$name 
  rownames(Full_peakMatrix_reads) <- Full_peakMatrix_names
  
  return(Full_peakMatrix_reads)
}
