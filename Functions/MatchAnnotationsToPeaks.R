# Script to identify overlaps between two GRanges, and transfer metadata from the subject to query.
# Useful for matching peak annotations e.g. motif co-ordinates with ArchR peaks and corresponding
# metadata such as nearby gene

MatchAnnotationsToPeaks <- function(query, subject) {
  
  Overlaps <- findOverlaps(query, subject)
  
  # Extract hits from overlaps
  query_hits <- queryHits(Overlaps)
  subject_hits <- subjectHits(Overlaps)
  
  # Extract the metadata from the subject GRanges object
  subject_metadata <- mcols(subject)
  
  # Subset the subject metadata based on subjectHits
  subject_metadata_subset <- subject_metadata[subject_hits, , drop = FALSE]
  
  # Subset the query GRanges object based on the query hits
  query_subset <- Sox8_motif_positions[query_hits]
  
  # Combine the metadata to the query GRanges object
  combined_metadata <- cbind(mcols(query_subset), subject_metadata_subset)
  
  # Assign the combined metadata back to the query GRanges object
  mcols(query_subset) <- combined_metadata
  
  return(query_subset)
}
