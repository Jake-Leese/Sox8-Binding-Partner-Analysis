Match_GRanges_to_Metadata <- function(query, subject) {
  
  subject_ranges_from_metadata <- GRanges(
    seqnames = subject$seqnames,   # Assuming 'seqnames' is a column in B's metadata
    ranges = IRanges(
      start = subject$start,       # Assuming 'start' is a column in B's metadata
      end = subject$end            # Assuming 'end' is a column in B's metadata
    )
  )
  
  overlaps <- findOverlaps(subject_ranges_from_metadata, query)
  
  query_subset <- query[subjectHits(overlaps)]
  
  return(query_subset)
}
