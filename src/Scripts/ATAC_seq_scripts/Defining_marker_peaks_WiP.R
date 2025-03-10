# Script to extract marker peaks for PPR and NC from full peakset. Output is NC, PPR, and Co-accessible peaksets. These are then used for downstream motif analysis and enhancer annotations

.libPaths("/R/libs/AT_ArcR_macs2")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
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

# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

# Load custom R functions
source("/data/Sox8_binding_partner_analysis/src/Functions/ArchRAddUniqueIdsToSe.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/ArchR_ExtractIds.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/AddArchRMetaData.R")


##############################################################################################################################
######################################### Get Marker Features against all ####################################################

ss4_markerPeaks <- getMarkerFeatures(
  ArchRProj = ss4ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

ss8_markerPeaks <- getMarkerFeatures(
  ArchRProj = ss8ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

# Plot heatmap of marker peaks by cell type
ss4_heatmapPeaks <- markerHeatmap(
  seMarker = ss4_markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE)
draw(ss4_heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

ss8_heatmapPeaks <- markerHeatmap(
  seMarker = ss8_markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE)
draw(ss8_heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# Extract GRanges list of marker peaks for all cell types. Using relatively relaxed thresholds: Log Fold change > 0, and FDR <= 0.1
ss4_markerList <- getMarkers(ss4_markerPeaks, cutOff = "Log2FC > 1 & FDR <= 0.1", returnGR = TRUE)
ss4_Placodal_markers <- ss4_markerList$Placodal
ss4_NC_markers <- ss4_markerList$NC

ss8_markerList <- getMarkers(ss8_markerPeaks, cutOff = "Log2FC > 1 & FDR <= 0.1", returnGR = TRUE)
ss8_Placodal_markers <- ss8_markerList$Placodal
ss8_NC_markers <- ss8_markerList$NC


##############################################################################################################################
##################################### Get Marker Features from pairwise tests ################################################


# Alternative method, rather than calculating Log2FC and FDR for each cell type against the all other cells as background, I will
# See if a series of pairwise tests between each combination of cell types produces more peaks for our cell types of interest.
# The rationale here, is that as long as a peak is comparatively more accessible than at least one other cell type, then it is likely to be open in the first cell type
# Similar method has been used before in this paper: https://www.nature.com/articles/nature22367#Sec9
# Will test this between broad ectodermal cell types

# Define the cell types present in your dataset
cell_types <- unique(ss4ArchRProj$transferred_scHelper_cell_type_broad)

# Output directory
output_dir <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/"
dir.create(output_dir, showWarnings = FALSE)

# Loop over each pair of cell types
for (i in seq_along(cell_types)) {
  for (j in seq_along(cell_types)) {
    if (i != j) {  # Avoid self-comparisons
      
      # Define groups
      cell_type_A <- cell_types[i]
      cell_type_B <- cell_types[j]
      
      # Print progress
      message("Comparing: ", cell_type_A, " vs ", cell_type_B)
      
      # Run getMarkerFeatures() for pairwise comparison
      markers <- getMarkerFeatures(
        ArchRProj = ss4ArchRProj,
        useMatrix = "PeakMatrix",
        groupBy = "transferred_scHelper_cell_type_broad",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon",
        useGroups = cell_type_A,  # Cell type of interest
        bgdGroups = cell_type_B   # Background cell type
      )
      
      # Extract peaks withLog2FC and FDR values
      markers_gr <- getMarkers(markers,
                 cutOff = "FDR <= 0.05 & Log2FC >= 1",
                 returnGR = TRUE)
        
      # Save results to file
      saveRDS(markers_gr, file = paste0(output_dir, cell_type_A, "_vs_", cell_type_B, "_accessible_peaks_05FDR.gr.rds"))
    }
  }
}

# Read in relevant GRanges

# ss4 PPR
#ss4_PPR_v_contam_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4Placodal_vs_Contam_accessible_peaks.gr.rds")$Placodal
#ss4_PPR_v_NC_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4Placodal_vs_NC_accessible_peaks.gr.rds")$Placodal
#ss4_PPR_v_Neural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4Placodal_vs_Neural_accessible_peaks.gr.rds")$Placodal
#ss4_PPR_v_NonNeural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4Placodal_vs_Non-neural_accessible_peaks.gr.rds")$Placodal   # Empty GRanges

# ss4 NC
#ss4_NC_v_contam_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4NC_vs_Contam_accessible_peaks.gr.rds")$NC
#ss4_NC_v_Placodal_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4NC_vs_Placodal_accessible_peaks.gr.rds")$NC
#ss4_NC_v_Neural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4NC_vs_Neural_accessible_peaks.gr.rds")$NC
#ss4_NC_v_NonNeural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/ss4NC_vs_Non-neural_accessible_peaks.gr.rds")$NC   # Empty GRanges

# ss8 PPR
#ss8_PPR_v_contam_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/ss8Placodal_vs_Contam_accessible_peaks.gr.rds")$Placodal
#ss8_PPR_v_NC_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/ss8Placodal_vs_NC_accessible_peaks.gr.rds")$Placodal
#ss8_PPR_v_Neural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/ss8Placodal_vs_Neural_accessible_peaks.gr.rds")$Placodal

# ss8 NC
#ss8_NC_v_contam_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/ss8NC_vs_Contam_accessible_peaks.gr.rds")$NC
#ss8_NC_v_Placodal_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/ss8NC_vs_Placodal_accessible_peaks.gr.rds")$NC
#ss8_NC_v_Neural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/ss8NC_vs_Neural_accessible_peaks.gr.rds")$NC

# Combine GRanges of differentially accessible peaks for each stage and cell type. C
#ss4_PPR_all_Peaks <- GenomicRanges::reduce(c(ss4_PPR_v_contam_GR, ss4_PPR_v_NC_GR, ss4_PPR_v_Neural_GR)) # The reduce() function combines GRanges while avoiding duplicates
#ss4_NC_all_Peaks <- GenomicRanges::reduce(c(ss4_NC_v_contam_GR, ss4_NC_v_Placodal_GR, ss4_NC_v_Neural_GR))

#ss8_PPR_all_Peaks <- GenomicRanges::reduce(c(ss8_PPR_v_contam_GR, ss8_PPR_v_NC_GR, ss8_PPR_v_Neural_GR)) 
#ss8_NC_all_Peaks <- GenomicRanges::reduce(c(ss8_NC_v_contam_GR, ss8_NC_v_Placodal_GR, ss8_NC_v_Neural_GR))

# Sort the GRanges by seqnames and ranges. sortSeqlevels deals with the character format, i.e. "chr*", and enables numerical sorting of seqnames 
#ss4_PPR_all_Peaks <- sort(sortSeqlevels(ss4_PPR_all_Peaks))
#ss4_NC_all_Peaks <- sort(sortSeqlevels(ss4_NC_all_Peaks))
#ss8_PPR_all_Peaks <- sort(sortSeqlevels(ss8_PPR_all_Peaks))
#ss8_NC_all_Peaks <- sort(sortSeqlevels(ss8_NC_all_Peaks))

# Performed these two commands with each combined GRanges object and confirmed that all peaks are accounted for and that there are no duplicates
all(ss4_NC_v_contam_GR %over% ss4_NC_all_Peaks) & all(ss4_NC_v_Placodal_GR %over% ss4_NC_all_Peaks) & all(ss4_NC_v_Neural_GR %over% ss4_NC_all_Peaks)  # Should return TRUE if all original ranges are accounted for
any(duplicated(paste0(seqnames(ss4_NC_all_Peaks), ":", ranges(ss4_NC_all_Peaks)))) # This checks for any duplicates (seqname and ranges) and returns a logical value


##############################################################################################################################
####################################### Extract ranges for downstream analysis ###############################################


# Check for overlap between PPR and NC only (based on pairwise comparisons)
all(ss4_PPR_v_NC_GR %over% ss4_NC_v_Placodal_GR) # FALSE
all(ss8_PPR_v_NC_GR %over% ss8_NC_v_Placodal_GR) # FALSE
# FALSE for both stages, therefore these can be used in downstream analysis for PPR/NC only (FDR < 0.1; Log2FC >= 1)

# Find overlapping ranges that are accessible in both cell types
ss4_PPR_open_non_differential <- ss4_PPR_all_Peaks[!(ss4_PPR_all_Peaks %over% ss4_PPR_v_NC_GR)] # Subsets out the ss4_PPR_all_Peaks that are differentially accessible compared to NC (i.e. found in ss4_PPR_NC_GR)
ss4_NC_open_non_differential <- ss4_NC_all_Peaks[!(ss4_NC_all_Peaks %over% ss4_NC_v_Placodal_GR)]

all(ss4_PPR_open_non_differential %over% ss4_PPR_v_NC_GR) # FALSE: No overlap
all(ss4_NC_open_non_differential %over% ss4_NC_v_Placodal_GR) # FALSE

sum(ss4_PPR_open_non_differential %over% ss4_NC_open_non_differential) # 2362 peaks are accessible across both cell types. Leaves the majority of "open" enhancers unique to each cell type, despite not being found in a pairwise differential accessibility analysis between the two
                                                                       # This is because with this method, some peaks that are differentially accessible between NC/PPR and other cell types, may not necessarily be differentially accessible between NC and PPR
sum(ss4_Placodal_markers %over% ss4_PPR_v_NC_GR)

#ss4_PPR_all_Peaks[ranges(ss4_PPR_all_Peaks) %in% ranges(ss4_PPR_all_Peaks)[duplicated(ranges(ss4_PPR_all_Peaks))]] # Command to see all duplicates in a GRanges, including the first


##############################################################################################################################
###################################### Comparing outputs from different methods ##############################################


# Comparing the different methods that we have used to get hold of "open" peaks in NC and PPR (based on FDR < 0.1)

# Using marker peaks to identify accessible peaks across both cell types
ss4_Placodal_markers_non_diff <- ss4_Placodal_markers[!(ss4_Placodal_markers %over% ss4_PPR_v_NC_GR)] # 4013 to 2406 peaks
ss4_NC_markers_non_diff <- ss4_NC_markers[!(ss4_NC_markers %over% ss4_NC_v_Placodal_GR)] # 616 to 259 peaks
sum(ss4_Placodal_markers_non_diff %over% ss4_NC_markers_non_diff) # Only 2 peaks with shared accessibility using this method

ss8_Placodal_markers_non_diff <- ss8_Placodal_markers[!(ss8_Placodal_markers %over% ss8_PPR_v_NC_GR)] # 13612 to 6212 peaks
ss8_NC_markers_non_diff <- ss8_NC_markers[!(ss8_NC_markers %over% ss8_NC_v_Placodal_GR)] # 4057 to 1657 peaks
sum(ss8_Placodal_markers_non_diff %over% ss8_NC_markers_non_diff) # Only 121 peaks with shared accessibility using this method

# Using outputs of individual pairwise DE analyses to identify accessible peaks
ss4_PPR_pairwise_non_diff <- ss4_PPR_all_Peaks[!(ss4_PPR_all_Peaks %over% ss4_PPR_v_NC_GR)] # 14276 to 9818 peaks
ss4_NC_pairwise_non_diff <- ss4_NC_all_Peaks[!(ss4_NC_all_Peaks %over% ss4_NC_v_Placodal_GR)] # 13148 to 8137 peaks
sum(ss4_PPR_pairwise_non_diff %over% ss4_NC_pairwise_non_diff) # 2362 peaks with shared accessibility

ss8_PPR_pairwise_non_diff <- ss8_PPR_all_Peaks[!(ss8_PPR_all_Peaks %over% ss8_PPR_v_NC_GR)] # 29685 to 15153 peaks
ss8_NC_pairwise_non_diff <- ss8_NC_all_Peaks[!(ss8_NC_all_Peaks %over% ss8_NC_v_Placodal_GR)] # 22396 to 11236 peaks
sum(ss8_PPR_pairwise_non_diff %over% ss8_NC_pairwise_non_diff) # 3342 peaks with shared accessibility

#################################################### CONCLUSION ##############################################################
# Defining peak accessibility based on differentially accessible between at least one cell-cell comparison gives a much higher number of peaks
# There is very little overlap between marker peak for NC and PPR, which means this is an inappropriate method for identifying shared accessibility



##############################################################################################################################
################################### Generating Final output files using FDR < 0.05 ###########################################

# Read in relevant GRanges

# ss4 PPR
ss4_PPR_v_contam_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/Placodal_vs_Contam_accessible_peaks_05FDR.gr.rds")$Placodal
ss4_PPR_v_NC_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/Placodal_vs_NC_accessible_peaks_05FDR.gr.rds")$Placodal
ss4_PPR_v_Neural_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/Placodal_vs_Neural_accessible_peaks_05FDR.gr.rds")$Placodal
#ss4_PPR_v_NonNeural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/Placodal_vs_Non-neural_accessible_peaks_05FDR.gr.rds")$Placodal   # Empty GRanges

# ss4 NC
ss4_NC_v_contam_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/NC_vs_Contam_accessible_peaks_05FDR.gr.rds")$NC
ss4_NC_v_Placodal_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/NC_vs_Placodal_accessible_peaks_05FDR.gr.rds")$NC
ss4_NC_v_Neural_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/NC_vs_Neural_accessible_peaks_05FDR.gr.rds")$NC
#ss4_NC_v_NonNeural_GR <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss4/NC_vs_Non-neural_accessible_peaks_05FDR.gr.rds")$NC   # Empty GRanges

# ss8 PPR
ss8_PPR_v_contam_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/Placodal_vs_Contam_accessible_peaks_05FDR.gr.rds")$Placodal
ss8_PPR_v_NC_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/Placodal_vs_NC_accessible_peaks_05FDR.gr.rds")$Placodal
ss8_PPR_v_Neural_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/Placodal_vs_Neural_accessible_peaks_05FDR.gr.rds")$Placodal

# ss8 NC
ss8_NC_v_contam_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/NC_vs_Contam_accessible_peaks_05FDR.gr.rds")$NC
ss8_NC_v_Placodal_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/NC_vs_Placodal_accessible_peaks_05FDR.gr.rds")$NC
ss8_NC_v_Neural_GR_0.05 <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/ss8/NC_vs_Neural_accessible_peaks_05FDR.gr.rds")$NC

# Check for overlap between PPR and NC only (based on pairwise comparisons)
all(ss4_PPR_v_NC_GR_0.05 %over% ss4_NC_v_Placodal_GR_0.05) # FALSE
all(ss8_PPR_v_NC_GR_0.05 %over% ss8_NC_v_Placodal_GR_0.05) # FALSE

# Combine GRanges of differentially accessible peaks for each stage and cell type. C
ss4_PPR_all_Peaks_0.05 <- GenomicRanges::reduce(c(ss4_PPR_v_contam_GR_0.05, ss4_PPR_v_NC_GR_0.05, ss4_PPR_v_Neural_GR_0.05)) # The reduce() function combines GRanges while avoiding duplicates
ss4_NC_all_Peaks_0.05 <- GenomicRanges::reduce(c(ss4_NC_v_contam_GR_0.05, ss4_NC_v_Placodal_GR_0.05, ss4_NC_v_Neural_GR_0.05))

ss8_PPR_all_Peaks_0.05 <- GenomicRanges::reduce(c(ss8_PPR_v_contam_GR_0.05, ss8_PPR_v_NC_GR_0.05, ss8_PPR_v_Neural_GR_0.05)) 
ss8_NC_all_Peaks_0.05 <- GenomicRanges::reduce(c(ss8_NC_v_contam_GR_0.05, ss8_NC_v_Placodal_GR_0.05, ss8_NC_v_Neural_GR_0.05))

# Sort the GRanges by seqnames and ranges. sortSeqlevels deals with the character format, i.e. "chr*", and enables numerical sorting of seqnames 
ss4_PPR_all_Peaks_0.05 <- sort(sortSeqlevels(ss4_PPR_all_Peaks_0.05))
ss4_NC_all_Peaks_0.05 <- sort(sortSeqlevels(ss4_NC_all_Peaks_0.05))
ss8_PPR_all_Peaks_0.05 <- sort(sortSeqlevels(ss8_PPR_all_Peaks_0.05))
ss8_NC_all_Peaks_0.05 <- sort(sortSeqlevels(ss8_NC_all_Peaks_0.05))

# Save peaksets
path <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/"
saveRDS(ss4_PPR_all_Peaks_0.05, paste0(path, "ss4_PPR_all_Peaks.gr.rds"))
saveRDS(ss4_NC_all_Peaks_0.05, paste0(path, "ss4_NC_all_Peaks.gr.rds"))
saveRDS(ss8_PPR_all_Peaks_0.05, paste0(path, "ss8_PPR_all_Peaks.gr.rds"))
saveRDS(ss8_NC_all_Peaks_0.05, paste0(path, "ss8_NC_all_Peaks.gr.rds"))

# Performed these two commands with each combined GRanges object and confirmed that all peaks are accounted for and that there are no duplicates
all(ss8_NC_v_contam_GR_0.05 %over% ss8_NC_all_Peaks_0.05) & all(ss8_NC_v_Placodal_GR_0.05 %over% ss8_NC_all_Peaks_0.05) & all(ss8_NC_v_Neural_GR_0.05 %over% ss8_NC_all_Peaks_0.05)  # Should return TRUE if all original ranges are accounted for
any(duplicated(paste0(seqnames(ss8_NC_all_Peaks_0.05), ":", ranges(ss8_NC_all_Peaks_0.05)))) # This checks for any duplicates (seqname and ranges) and returns a logical value

# Extracting peaks with shared accessibility for NC and PPR
ss4_PPR_pairwise_non_diff_0.05 <- ss4_PPR_all_Peaks_0.05[!(ss4_PPR_all_Peaks_0.05 %over% ss4_PPR_v_NC_GR_0.05)] # 9833 to 6955 peaks
ss4_NC_pairwise_non_diff_0.05 <- ss4_NC_all_Peaks_0.05[!(ss4_NC_all_Peaks_0.05 %over% ss4_NC_v_Placodal_GR_0.05)] # 8932 to 5596 peaks
sum(ss4_PPR_pairwise_non_diff_0.05 %over% ss4_NC_pairwise_non_diff_0.05) # 1576 peaks with shared accessibility
ss4_PPR_NC_Shared <- ss4_PPR_pairwise_non_diff_0.05[ss4_PPR_pairwise_non_diff_0.05 %over% ss4_NC_pairwise_non_diff_0.05]

ss8_PPR_pairwise_non_diff_0.05 <- ss8_PPR_all_Peaks_0.05[!(ss8_PPR_all_Peaks_0.05 %over% ss8_PPR_v_NC_GR_0.05)] # 22347 to 11117 peaks
ss8_NC_pairwise_non_diff_0.05 <- ss8_NC_all_Peaks_0.05[!(ss8_NC_all_Peaks_0.05 %over% ss8_NC_v_Placodal_GR_0.05)] # 16117 to 8304 peaks
sum(ss8_PPR_pairwise_non_diff_0.05 %over% ss8_NC_pairwise_non_diff_0.05) # 2318 peaks with shared accessibility
ss8_PPR_NC_Shared <- ss8_PPR_pairwise_non_diff_0.05[ss8_PPR_pairwise_non_diff_0.05 %over% ss8_NC_pairwise_non_diff_0.05]

path <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/"
saveRDS(ss4_PPR_v_NC_GR_0.05, paste0(path, "ss4_PPR_only.gr.rds"))
saveRDS(ss4_NC_v_Placodal_GR_0.05, paste0(path, "ss4_NC_only.gr.rds"))
saveRDS(ss4_PPR_NC_Shared, paste0(path, "ss4_PPR_and_NC_shared.gr.rds"))
saveRDS(ss8_PPR_v_NC_GR_0.05, paste0(path, "ss8_PPR_only.gr.rds"))
saveRDS(ss8_NC_v_Placodal_GR_0.05, paste0(path, "ss8_NC_only.gr.rds"))
saveRDS(ss8_PPR_NC_Shared, paste0(path, "ss8_PPR_and_NC_shared.gr.rds"))

