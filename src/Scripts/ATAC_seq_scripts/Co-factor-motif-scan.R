# This script extracts the Genomic positions of a transcription factor binding motif (in this case Sox8), and using the R package, motifmatchR,
# to search within a specified range of these positions for any neighbouring transcription factor binding motifs.

.libPaths("/scratch/users/k2370243/software/R/4.2")
setwd('/scratch/users/k2370243/scATACseq/')

# Dependencies for this script
library(ArchR)
library(JASPAR2020)
library(BSgenome.Ggallus.UCSC.galGal6)
library(motifmatchr)
library(GenomicRanges)
library(TFBSTools)
library(dplyr)
library(tidyr)
library(purrr)
library(remotes)

# BiocManager:

# Loading ArchR project
ss8_Sox8_ArchRProj <- loadArchRProject(path = '/scratch/users/k2370243/scATACseq/ss8_Sox8_Save_ArchR/', force = FALSE)

# Extracting the Sox8 motif positions from the subsetted ArchR project.
# This returns a GRanges List object
Sox8_motifPositions <- getPositions(ss8_Sox8_ArchRProj, annoName = "SOX8")

# Convert the GRanges List object to a GRanges object
Sox8_PositionsGR <- as(Sox8_motifPositions$SOX8, "GRanges")

# Subset Sox8 Granges by placodal 
Placodal_overlaps <- findOverlaps(Sox8_PositionsGR, Placodal_peaks)

Placodal_Sox8_Positions <- Sox8_PositionsGR[queryHits(Placodal_overlaps)]

# Subset Sox8 Granges by NC 
NC_overlaps <- findOverlaps(Sox8_PositionsGR, NC_Peaks)

NC_Sox8_Positions <- Sox8_PositionsGR[queryHits(NC_overlaps)]

# Function to expand the ranges of a GRanges object
expand_ranges_upstream <- function(gr, expansion_amount) {
  
  # Extract the GRanges object as a DataFrame
  ranges_df <- as.data.frame(gr)
  
  # Expand the ranges (specified by start and end columns) by the specified amount
  ranges_df$start <- ranges_df$start - expansion_amount
  ranges_df$end <- ranges_df$start
  
  # Create a new GRanges object with expanded ranges. Keep extra columns as metadata.
  expanded_gr <- makeGRangesFromDataFrame(ranges_df, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
  
  return(expanded_gr)
}

# Function to expand the ranges of a GRanges object
expand_ranges_downstream <- function(gr, expansion_amount) {
  
  # Extract the GRanges object as a DataFrame
  ranges_df <- as.data.frame(gr)
  
  # Expand the ranges (specified by start and end columns) by the specified amount
  ranges_df$start <- ranges_df$end
  ranges_df$end <- ranges_df$end + expansion_amount
  
  # Create a new GRanges object with expanded ranges. Keep extra columns as metadata.
  expanded_gr <- makeGRangesFromDataFrame(ranges_df, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
  
  return(expanded_gr)
}


# Create a new GRanges object using expand_ranges() with the specified range expansion as the second argument
Sox8_Placodal_upstream_50bp <- expand_ranges_upstream(Placodal_Sox8_Positions, 50)
Sox8_Placodal_downstream_50bp <- expand_ranges_downstream(Placodal_Sox8_Positions, 50)

Sox8_Placodes_50bp_expanded_excl_Sox8motif <- c(Sox8_Placodes_upstream_50bp, Sox8_Placodes_downstream_50bp)

# Create a new GRanges object using expand_ranges() with the specified range expansion as the second argument
Sox8_NC_upstream_50bp <- expand_ranges_upstream(NC_Sox8_Positions, 50)
Sox8_NC_downstream_50bp <- expand_ranges_downstream(NC_Sox8_Positions, 50)

Sox8_NC_50bp_expanded_excl_Sox8motif <- c(Sox8_NC_upstream_50bp, Sox8_NC_downstream_50bp)

# download motif database
motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# rename each motif to have TF name
name_vector <- c()
for (i in 1:length(motifList)){
  name <- name(motifList[[i]])
  name_vector <- c(name_vector, name)
}
names(motifList) <- name_vector

# Find the nucleotide frequencies of the original peakset in order to replicate the ArchR addMotifAnnotation() analysis
# Get peakset from ss8_ArchR project. 
Full_ss8_peakSet <- getPeakSet(ArchRProj = ss8ArchRProj)
# Get the ACGT frequencies from the full ss8 peakset
freqs <- alphabetFrequency(getSeq(BSgenome.Ggallus.UCSC.galGal6, Full_ss8_peakSet))
ACGT_freqs <- c(sum(freqs[,"A"])/sum(freqs),
                sum(freqs[,"C"])/sum(freqs),
                sum(freqs[,"G"])/sum(freqs),
                sum(freqs[,"T"])/sum(freqs))
# ss8_peakset bg ACTG_freqs = 0.2487819, 0.2512324, 0.2512223, 0.2487634

# Scan GRanges object for motifs
NC_motif_matches_scores <- matchMotifs(
  motifList, 
  Sox8_NC_50bp_expanded_excl_Sox8motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

# Scan GRanges object for motifs and return motif positions
motif_matches_positions <- matchMotifs(
  motifList, 
  Sox8_20bp_expanded_excl_Sox8motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "positions", 
  p.cutoff = 1e-05)

# Extract iundividual motif scores and motif matches from the summarizes experiment
motif_scores <- assays(motif_matches_scores)$motifScores
motif_matches <- assays(NC_motif_matches_scores)$motifMatches

# Calculate the sum of TRUE values for each collumn of the motif match matrix
column_sums <- colSums(motif_matches)

# Order the transcription factors based on the sums in descending order
ordered_columns <- names(sort(column_sums, decreasing = TRUE))

# Use the ordered transcription factors, to sort the column_sums named vector
NC_sorted_TF_matches <- column_sums[ordered_columns]

# Display the ordered columns
print(ordered_columns)

SOX8_motifpositions_new <- motif_matches_positions$SOX8

countOverlaps(Sox8_PositionsGR,
            SOX8_motifpositions_new)
