# Function to take add peak matrix metadata to corresponding peaks in a GRanges object

AddArchRMetaData <- function(GRanges_object, ArchR_project) {

    tmp_peaks = data.frame(ArchR_project@peakSet) # Make a dataframe from the ArchR project peakset
    tmp_diff_peaks = data.frame(GRanges_object)   # Make a dataframe from the GRanges object
    diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, # Join the dataframes based on seqnames and Iranges, maintaining order of tmp_diff_peaks
                                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
    diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_") # Add a new column 
    mcols(GRanges_object) = diff_peaks_join_peakset

  return(GRanges_object)
}
