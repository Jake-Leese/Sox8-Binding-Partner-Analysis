# Function to expand the ranges a GRanges object by a specific value
# The function adds the specified range upstream and downstream of the outer values of the existing range
# Importantly, it removes the original range, meaning that downstream motifscanning only focuses on the surrounding sequence
# As a result, this doubles the number of ranges in our GRanges object
# Function to expand the ranges of a GRanges object
expand_ranges <- function(GRanges, expansion_amount, KeepExistingRange) {
  
  if (KeepExistingRange == TRUE) {
    ranges <- as.data.frame(GRanges)
    ranges$start <- ranges$start - expansion_amount
    ranges$end <- ranges$end + expansion_amount
    expanded_gr <- makeGRangesFromDataFrame(ranges, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
    return(expanded_gr)
  }
  
  if (KeepExistingRange == FALSE) {
    
    ranges_upstream <- as.data.frame(GRanges)
    ranges_upstream$end <- ranges_upstream$start
    ranges_upstream$start <- ranges_upstream$start - expansion_amount
    ranges_upstream <- makeGRangesFromDataFrame(ranges_upstream, 
                                                keep.extra.columns = TRUE, 
                                                start.field = "start", 
                                                end.field = "end")
    
    ranges_downstream <- as.data.frame(GRanges)
    ranges_downstream$start <- ranges_downstream$end
    ranges_downstream$end <- ranges_downstream$end + expansion_amount
    ranges_downstream <- makeGRangesFromDataFrame(ranges_downstream, 
                                                keep.extra.columns = TRUE, 
                                                start.field = "start", 
                                                end.field = "end")
    
    Expanded_GR <- c(ranges_upstream, ranges_downstream)
    
    return(Expanded_GR)
  }
}


