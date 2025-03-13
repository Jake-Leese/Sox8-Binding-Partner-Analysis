# Script to:
# 1. extract PPR and NC accessible peaks from full peakset. 
# 2. Scan peaks for Sox8 and subset accordingly. 
# 3. Plot differential accessibility between PPR and NC Sox8 peaks
# 4. Output Sox8+ differentially accessible and mutually accessible peaks for NC and SOX8 at ss4 and ss8 that are ready for SpaMo analysis
# Ran using alexthiery-schelper-archr_dev_macs2-schelper-0.3.5.img containers

.libPaths("/R/libs/AT_ArcR_macs2")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(TFBSTools)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(gtools)
library(parallel)
library(clustree)
library(ComplexHeatmap)
library(BSgenome.Ggallus.UCSC.galGal6)
library(scHelper)
library(ggrepel)
library(BSgenome.Ggallus.UCSC.galGal6)
library(TxDb.Ggallus.UCSC.galGal6.refGene)
library(org.Gg.eg.db)
library(JASPAR2020)
library(motifmatchr)
library(universalmotif)

# Load custom R functions
source("/data/Sox8_binding_partner_analysis/src/Functions/MatchAnnotationsToPeaks.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/expand_ranges.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/ArchRAddUniqueIdsToSe.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/ArchR_ExtractIds.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/AddArchRMetaData.R")

# Set data paths
path_GR <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/March_2025/"
path_MEME <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/For_MEME/"

# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

##############################################################################################################################
############################ Define NC and PPR accessible peaks using pairwise tests ########################################

# Function to perform pairwise tests between defined cell types. Returns a GRanges list of all the pairwise comparisons

Extract_pairwise_marker_peaks <- function(ArchRProj, cell_type, meta_column) {
  
  # Get unique cell types from the specified metadata column
  cell_types <- unique(ArchRProj@cellColData[[meta_column]])
  
  # Check if the specified Cell_type exists in the metadata
  if (!(cell_type %in% cell_types)) {
    stop("Specified cell_type not found in the metadata column.")
  }
  
  # Remove the selected cell_type from the list to perform pairwise comparisons
  other_cell_types <- setdiff(cell_types, cell_type)
  
  # Store results in a list
  marker_results <- GRangesList()
  
  for (cell_type_B in other_cell_types) {
    message(paste("Performing differential accessibility analysis:", cell_type, "vs", cell_type_B))
    
    # Define groups for pairwise comparison
    groups <- list(cell_type = cell_type, Other = cell_type_B)
    
    # Perform differential accessibility analysis
    markers <- getMarkerFeatures(
      ArchRProj = ArchRProj,
      useMatrix = "PeakMatrix",
      groupBy = meta_column,
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      useGroups = cell_type,
      bgdGroups = cell_type_B
    )
    
    # Extract significant peaks (FDR <= 0.05, Log2FC >= 1)
    marker_list <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
    
    # Convert to GRanges and store in GRangesList
    if (!is.null(marker_list[[cell_type]])) {
      marker_results[[paste0(cell_type, "_vs_", cell_type_B)]] <- marker_list[[cell_type]]
    }
  }
  return(marker_results)
}



ss4_PPR_accessible_peaks <- Extract_pairwise_marker_peaks(ss4ArchRProj, "Placodal", "transferred_scHelper_cell_type_broad")
ss4_NC_accessible_peaks <- Extract_pairwise_marker_peaks(ss4ArchRProj, "NC", "transferred_scHelper_cell_type_broad")
ss8_PPR_accessible_peaks <- Extract_pairwise_marker_peaks(ss8ArchRProj, "Placodal", "transferred_scHelper_cell_type_broad")
ss8_NC_accessible_peaks <- Extract_pairwise_marker_peaks(ss8ArchRProj, "NC", "transferred_scHelper_cell_type_broad")

# Output directory
output_dir <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/"

# Save GRanges_lists
saveRDS(ss4_PPR_accessible_peaks, paste0(output_dir, "ss4/PPR_pairwise_peaks.gr.list.rds"))
saveRDS(ss4_NC_accessible_peaks, paste0(output_dir, "ss4/NC_pairwise_peaks.gr.list.rds"))
saveRDS(ss8_PPR_accessible_peaks, paste0(output_dir, "ss8/PPR_pairwise_peaks.gr.list.rds"))
saveRDS(ss8_NC_accessible_peaks, paste0(output_dir, "ss8/NC_pairwise_peaks.gr.list.rds"))

# Read back in GRanges Lists
ss4_PPR_accessible_peaks <- readRDS(paste0(output_dir, "ss4/PPR_pairwise_peaks.gr.list.rds"))
ss4_NC_accessible_peaks <- readRDS(paste0(output_dir, "ss4/NC_pairwise_peaks.gr.list.rds"))
ss8_PPR_accessible_peaks <- readRDS(paste0(output_dir, "ss8/PPR_pairwise_peaks.gr.list.rds"))
ss8_NC_accessible_peaks <- readRDS(paste0(output_dir, "ss8/NC_pairwise_peaks.gr.list.rds"))


# Consolidate all ss4 and ss8 PPR and NC accessible peaks
# Function to unify GRanges within a GRanges list and avoid duplicates
Unify_GR_list <- function(GRanges_list) {
  combined_gr <- unlist(GRanges_list, use.names = FALSE) # Flattens GRangesList into single Granges object
  combined_gr <- unique(combined_gr) # Removes any duplicates
  combined_gr <- sort(sortSeqlevels(combined_gr)) # Sorting the GRanges object by ranges
  return(combined_gr)
}

ss4_PPR_accessible_peaks_combined <- Unify_GR_list(ss4_PPR_accessible_peaks)
ss4_NC_accessible_peaks_combined <- Unify_GR_list(ss4_NC_accessible_peaks)
ss8_PPR_accessible_peaks_combined <- Unify_GR_list(ss8_PPR_accessible_peaks)
ss8_NC_accessible_peaks_combined <- Unify_GR_list(ss8_NC_accessible_peaks)

# Combine all NC and PPR peaks from ss4 and ss8 for downstream SOX8 motif scanning
Combined_NC_PPR_Peaks <- sort(sortSeqlevels(unique(c(ss4_PPR_accessible_peaks_combined, 
                                                     ss4_NC_accessible_peaks_combined,
                                                     ss8_PPR_accessible_peaks_combined,
                                                     ss8_NC_accessible_peaks_combined))))

# Check for duplicated ranges
ss8_NC_accessible_peaks_combined[ranges(ss8_NC_accessible_peaks_combined) %in% 
                                  ranges(ss8_NC_accessible_peaks_combined)[duplicated(ranges(ss8_NC_accessible_peaks_combined))]]

# Extract peaks which are differentially accessible between PPR and NC specifically
ss4_PPR_only_peaks <- ss4_PPR_accessible_peaks$Placodal_vs_NC  # 2879
ss4_NC_only_peaks <- ss4_NC_accessible_peaks$NC_vs_Placodal    # 3340
ss8_PPR_only_peaks <- ss8_PPR_accessible_peaks$Placodal_vs_NC  # 11238
ss8_NC_only_peaks <- ss8_NC_accessible_peaks$NC_vs_Placodal    # 7821

# Ensure there is no overlap between PPR_only and NC_only for each stage
any(ss4_PPR_only_peaks %over% ss4_NC_only_peaks) # FALSE
any(ss8_PPR_only_peaks %over% ss8_NC_only_peaks) # FALSE

# Find mutually accessible peaks for NC and PPR at both stages
ss4_PPR_and_NC_peaks <- ss4_PPR_accessible_peaks_combined[ss4_PPR_accessible_peaks_combined %over% ss4_NC_accessible_peaks_combined]  # 1710 peaks
ss8_PPR_and_NC_peaks <- ss8_PPR_accessible_peaks_combined[ss8_PPR_accessible_peaks_combined %over% ss8_NC_accessible_peaks_combined]  # 2648 peaks

#########################################################################################################
######################################## MOTIF SCANNING #################################################


# Define background nucleotide frequencies (based on full ArchR peakset)
ACTG_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Get the SOX8 motif from JASPAR2020
SOX8_tfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM", name = "SOX8"))[[1]]

# Scan for motif occurrences in the peakset
Sox8_motif_hits <- matchMotifs(SOX8_tfm, Combined_NC_PPR_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                          out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
saveRDS(Sox8_motif_hits, paste0(path_GR, "Sox8_positions_all_peaks.gr.rds"))

Sox8_motif_hits <- readRDS(paste0(path_GR, "Sox8_positions_all_peaks.gr.rds"))

# Subset original GRanges object to only those with a SOX8 binding motif
NC_PPR_combined_SOX8_gr <- subsetByOverlaps(Combined_NC_PPR_Peaks, Sox8_motif_hits)

# Subsetting GRanges by ranges that have a Sox8 motif
# Making peak sets for:

# 1. Differentially accessible NC and PPR between each other
ss4_PPR_v_NC_Sox8 <- ss4_PPR_only_peaks[ss4_PPR_only_peaks %over% NC_PPR_combined_SOX8_gr]       # 296 peaks
ss4_NC_v_PPR_Sox8 <- ss4_NC_only_peaks[ss4_NC_only_peaks %over% NC_PPR_combined_SOX8_gr]         # 674 peaks
ss8_PPR_v_NC_Sox8 <- ss8_PPR_only_peaks[ss8_PPR_only_peaks %over% NC_PPR_combined_SOX8_gr]       # 1234 peaks
ss8_NC_v_PPR_Sox8 <- ss8_NC_only_peaks[ss8_NC_only_peaks %over% NC_PPR_combined_SOX8_gr]         # 1591 peaks

# 2. Mutually accessible peaks
ss4_NC_and_PPR_Sox8 <- ss4_PPR_and_NC_peaks[ss4_PPR_and_NC_peaks %over% NC_PPR_combined_SOX8_gr] # 338 peaks
ss8_NC_and_PPR_Sox8 <- ss8_PPR_and_NC_peaks[ss8_PPR_and_NC_peaks %over% NC_PPR_combined_SOX8_gr] # 484 peaks

# 3. All PPR and NC accessible peaks (differentially accessible in at least one pairwise cell_type comparison)
ss4_PPR_all_Sox8 <- ss4_PPR_accessible_peaks_combined[ss4_PPR_accessible_peaks_combined %over% NC_PPR_combined_SOX8_gr]  # 1240 peaks
ss4_NC_all_Sox8 <- ss4_NC_accessible_peaks_combined[ss4_NC_accessible_peaks_combined %over% NC_PPR_combined_SOX8_gr]     # 1760 peaks
ss8_PPR_all_Sox8 <- ss8_PPR_accessible_peaks_combined[ss8_PPR_accessible_peaks_combined %over% NC_PPR_combined_SOX8_gr]  # 2751 peaks
ss8_NC_all_Sox8 <- ss8_NC_accessible_peaks_combined[ss8_NC_accessible_peaks_combined %over% NC_PPR_combined_SOX8_gr]     # 3114 peaks

# Save all Sox8 peak files

saveRDS(ss4_PPR_v_NC_Sox8, paste0(path_GR, "ss4_PPR_v_NC_Sox8.gr.rds"))
saveRDS(ss4_NC_v_PPR_Sox8, paste0(path_GR, "ss4_NC_v_PPR_Sox8.gr.rds"))
saveRDS(ss8_PPR_v_NC_Sox8, paste0(path_GR, "ss8_PPR_v_NC_Sox8.gr.rds"))
saveRDS(ss8_NC_v_PPR_Sox8, paste0(path_GR, "ss8_NC_v_PPR_Sox8.gr.rds"))
saveRDS(ss4_NC_and_PPR_Sox8, paste0(path_GR, "ss4_NC_and_PPR_Sox8.gr.rds"))
saveRDS(ss8_NC_and_PPR_Sox8, paste0(path_GR, "ss8_NC_and_PPR_Sox8.gr.rds"))
saveRDS(ss4_PPR_all_Sox8, paste0(path_GR, "ss4_PPR_all_Sox8.gr.rds"))
saveRDS(ss4_NC_all_Sox8, paste0(path_GR, "ss4_NC_all_Sox8.gr.rds"))
saveRDS(ss8_PPR_all_Sox8, paste0(path_GR, "ss8_PPR_all_Sox8.gr.rds"))
saveRDS(ss8_NC_all_Sox8, paste0(path_GR, "ss8_NC_all_Sox8.gr.rds"))


#########################################################################################################
################################### Preparing MEME outputs ##############################################

# Read in Sox8 peak files
lapply(list.files(path_GR, pattern = "\\.gr\\.rds$", full.names = TRUE), 
       function(f) assign(gsub("\\.gr\\.rds$", "", basename(f)), readRDS(f), envir = .GlobalEnv))

# load in all Sox8 motif co-ordinates
Sox8_motif_hits <- readRDS(paste0(path_GR, "Sox8_positions_all_peaks.gr.rds"))[[1]]

# Function to subset Sox8 binding co-ordinated (query) based on peakset (subset) of interest. Adds the peakset ranges as metadata for reference down the line
MatchSox8AnnotationsToPeaks <- function(query, subject) {

if(any(duplicated(subject)) == TRUE) {
  paste("Duplicates in subject GRanges")
}

Overlaps <- findOverlaps(query, subject)

# Extract hits from overlaps
query_hits <- queryHits(Overlaps)
subject_hits <- subjectHits(Overlaps)

# Extract the metadata from the subject GRanges object
subject_ranges <- DataFrame(seqnames = seqnames(subject), start = start(subject), end = end(subject))

# Subset/rearrange the subject metadata based on subjectHits
subject_ranges_subset <- subject_ranges[subject_hits, , drop = FALSE]

# Subset the query GRanges object based on the query hits
query_subset <- query[query_hits]

# Combine the metadata to the query GRanges object
combined_metadata <- cbind(mcols(query_subset), subject_ranges_subset)

# Assign the combined metadata back to the query GRanges object
mcols(query_subset) <- combined_metadata

return(query_subset)

}

# Subset Sox8 motif positions by Sox8 PPR/NC peak files. # Metadata columns are score (Sox8 motif match), seqnames, start and end (original ArchR peak ranges)
ss4_PPR_v_NC_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss4_PPR_v_NC_Sox8)
ss4_NC_v_PPR_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss4_NC_v_PPR_Sox8)
ss8_PPR_v_NC_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss8_PPR_v_NC_Sox8)
ss8_NC_v_PPR_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss8_NC_v_PPR_Sox8)
ss4_NC_and_PPR_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss4_NC_and_PPR_Sox8)
ss8_NC_and_PPR_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss8_NC_and_PPR_Sox8)
ss4_PPR_all_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss4_PPR_all_Sox8)
ss4_NC_all_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss4_NC_all_Sox8)
ss8_PPR_all_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss8_PPR_all_Sox8)
ss8_NC_all_Sox8_positions <- MatchSox8AnnotationsToPeaks(Sox8_motif_hits, ss8_NC_all_Sox8)


# Expand GRanges (+/- 100 bp) around Sox8 motif positions
ss4_PPR_v_NC_Sox8_positions_100bp_surr <- expand_ranges(ss4_PPR_v_NC_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_v_PPR_Sox8_positions_100bp_surr <- expand_ranges(ss4_NC_v_PPR_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_PPR_v_NC_Sox8_positions_100bp_surr <- expand_ranges(ss8_PPR_v_NC_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_v_PPR_Sox8_positions_100bp_surr <- expand_ranges(ss8_NC_v_PPR_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_and_PPR_Sox8_positions_100bp_surr <- expand_ranges(ss4_NC_and_PPR_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_and_PPR_Sox8_positions_100bp_surr <- expand_ranges(ss8_NC_and_PPR_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_PPR_all_Sox8_positions_100bp_surr <- expand_ranges(ss4_PPR_all_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_all_Sox8_positions_100bp_surr <- expand_ranges(ss4_NC_all_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_PPR_all_Sox8_positions_100bp_surr <- expand_ranges(ss8_PPR_all_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_all_Sox8_positions_100bp_surr <- expand_ranges(ss8_NC_all_Sox8_positions, expansion_amount = 100, KeepExistingRange = TRUE)


# Save as BED files for enhancer annotation

# Convert to a data frame
GRanges_to_Bed <- function(GRanges) {
  BED_df <- data.frame(
    chrom = gsub("^chr", "", as.character(seqnames(GRanges))),
    start = start(GRanges) - 1, # Convert to 0-based start for BED format
    end = end(GRanges),
    strand = as.character(strand(GRanges))
  )
  return(BED_df)
}
ss4_PPR_all_Sox8_positions_100bp_bed <- GRanges_to_Bed(ss4_PPR_all_Sox8_positions_100bp_surr)
ss4_NC_all_Sox8_positions_100bp_bed <- GRanges_to_Bed(ss4_NC_all_Sox8_positions_100bp_surr)
ss8_PPR_all_Sox8_positions_100bp_bed <- GRanges_to_Bed(ss8_PPR_all_Sox8_positions_100bp_surr)
ss8_NC_all_Sox8_positions_100bp_bed <- GRanges_to_Bed(ss8_NC_all_Sox8_positions_100bp_surr)

write.table(ss4_PPR_all_Sox8_positions_100bp_bed, file = "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/ss4_PPR_all_Sox8_positions_100bp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ss4_NC_all_Sox8_positions_100bp_bed, file = "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/ss4_NC_all_Sox8_positions_100bp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ss8_PPR_all_Sox8_positions_100bp_bed, file = "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/ss8_PPR_all_Sox8_positions_100bp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ss8_NC_all_Sox8_positions_100bp_bed, file = "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/ss8_NC_all_Sox8_positions_100bp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#### CREATE AND SAVE FASTA FILES FOR PEAKSETS AS INPUTS FOR MEME ####

GalGal6 <- BSgenome.Ggallus.UCSC.galGal6

# Get raw sequence using the GRanges objects
ss4_PPR_v_NC_Sox8_seq <- getSeq(GalGal6, ss4_PPR_v_NC_Sox8_positions_100bp_surr)
ss4_NC_v_PPR_Sox8_seq <- getSeq(GalGal6, ss4_NC_v_PPR_Sox8_positions_100bp_surr)
ss8_PPR_v_NC_Sox8_seq <- getSeq(GalGal6, ss8_PPR_v_NC_Sox8_positions_100bp_surr)
ss8_NC_v_PPR_Sox8_seq <- getSeq(GalGal6, ss8_NC_v_PPR_Sox8_positions_100bp_surr)
ss4_NC_and_PPR_Sox8_seq <- getSeq(GalGal6, ss4_NC_and_PPR_Sox8_positions_100bp_surr)
ss8_NC_and_PPR_Sox8_seq <- getSeq(GalGal6, ss8_NC_and_PPR_Sox8_positions_100bp_surr)
ss4_PPR_all_Sox8_seq <- getSeq(GalGal6, ss4_PPR_all_Sox8_positions_100bp_surr)
ss4_NC_all_Sox8_Sox8_seq <- getSeq(GalGal6, ss4_NC_all_Sox8_positions_100bp_surr)
ss8_PPR_all_Sox8_seq <- getSeq(GalGal6, ss8_PPR_all_Sox8_positions_100bp_surr)
ss8_NC_all_Sox8_Sox8_seq <- getSeq(GalGal6, ss8_NC_all_Sox8_positions_100bp_surr)

# Adding seqnames and ranges as unique sequence identifiers
names(ss4_PPR_v_NC_Sox8_seq) <- paste0(ss4_PPR_v_NC_Sox8_positions_100bp_surr@seqnames, ":", ss4_PPR_v_NC_Sox8_positions_100bp_surr@ranges)
names(ss4_NC_v_PPR_Sox8_seq) <- paste0(ss4_NC_v_PPR_Sox8_positions_100bp_surr@seqnames, ":", ss4_NC_v_PPR_Sox8_positions_100bp_surr@ranges)
names(ss8_PPR_v_NC_Sox8_seq) <- paste0(ss8_PPR_v_NC_Sox8_positions_100bp_surr@seqnames, ":", ss8_PPR_v_NC_Sox8_positions_100bp_surr@ranges)
names(ss8_NC_v_PPR_Sox8_seq) <- paste0(ss8_NC_v_PPR_Sox8_positions_100bp_surr@seqnames, ":", ss8_NC_v_PPR_Sox8_positions_100bp_surr@ranges)

names(ss4_NC_and_PPR_Sox8_seq) <- paste0(ss4_NC_and_PPR_Sox8_positions_100bp_surr@seqnames, ":", ss4_NC_and_PPR_Sox8_positions_100bp_surr@ranges)
names(ss8_NC_and_PPR_Sox8_seq) <- paste0(ss8_NC_and_PPR_Sox8_positions_100bp_surr@seqnames, ":", ss8_NC_and_PPR_Sox8_positions_100bp_surr@ranges)

names(ss4_PPR_all_Sox8_seq) <- paste0(ss4_PPR_all_Sox8_positions_100bp_surr@seqnames, ":", ss4_PPR_all_Sox8_positions_100bp_surr@ranges)
names(ss4_NC_all_Sox8_Sox8_seq) <- paste0(ss4_NC_all_Sox8_positions_100bp_surr@seqnames, ":", ss4_NC_all_Sox8_positions_100bp_surr@ranges)
names(ss8_PPR_all_Sox8_seq) <- paste0(ss8_PPR_all_Sox8_positions_100bp_surr@seqnames, ":", ss8_PPR_all_Sox8_positions_100bp_surr@ranges)
names(ss8_NC_all_Sox8_Sox8_seq) <- paste0(ss8_NC_all_Sox8_positions_100bp_surr@seqnames, ":", ss8_NC_all_Sox8_positions_100bp_surr@ranges)

# Write out sequences as .fasta files
writeXStringSet(ss4_PPR_v_NC_Sox8_seq, file=paste0(path_MEME, "ss4_PPR_v_NC_Sox8_seq.fasta"))
writeXStringSet(ss4_NC_v_PPR_Sox8_seq, file=paste0(path_MEME, "ss4_NC_v_PPR_Sox8_seq.fasta"))
writeXStringSet(ss8_PPR_v_NC_Sox8_seq, file=paste0(path_MEME, "ss8_PPR_v_NC_Sox8_seq.fasta"))
writeXStringSet(ss8_NC_v_PPR_Sox8_seq, file=paste0(path_MEME, "ss8_NC_v_PPR_Sox8_seq.fasta"))

writeXStringSet(ss4_NC_and_PPR_Sox8_seq, file=paste0(path_MEME, "ss4_NC_and_PPR_Sox8_seq.fasta"))
writeXStringSet(ss8_NC_and_PPR_Sox8_seq, file=paste0(path_MEME, "ss8_NC_and_PPR_Sox8_seq.fasta"))

writeXStringSet(ss4_PPR_all_Sox8_seq, file=paste0(path_MEME, "ss4_PPR_all_Sox8_seq.fasta"))
writeXStringSet(ss4_NC_all_Sox8_Sox8_seq, file=paste0(path_MEME, "ss4_NC_all_Sox8_Sox8_seq.fasta"))
writeXStringSet(ss8_PPR_all_Sox8_seq, file=paste0(path_MEME, "ss8_PPR_all_Sox8_seq.fasta"))
writeXStringSet(ss8_NC_all_Sox8_Sox8_seq, file=paste0(path_MEME, "ss8_NC_all_Sox8_Sox8_seq.fasta"))








markers <- getMarkerFeatures(
  ArchRProj = ss4ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "NC"
)

marker_list <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
marker_list[["Placodal"]]
