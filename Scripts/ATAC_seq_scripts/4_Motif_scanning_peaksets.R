# Script to:
# A. Prepare GRanges objects for each peakset that specify the sequence (of chosen length) surrounding the Sox8 binding motifs
# B. Prepare motif binding matrices for TFs of interest (based on scRNA-seq analysis)
# C. Scan Sox8-surrounding sequences for TFs of interest. Save the outputs (Motif Matches & Motif Positions) for downstream filtering and analysis

.libPaths("/R/libs/ArchR_Seurat_R_441")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(ArchR)
library(igraph)
library(Seurat)
library(TFBSTools)
library(dplyr)
library(tidyr)
library(purrr)
library(BSgenome.Ggallus.UCSC.galGal6)
library(TxDb.Ggallus.UCSC.galGal6.refGene)
library(org.Gg.eg.db)
library(JASPAR2024)
#### 

source("/data/Sox8_binding_partner_analysis/Functions/MatchAnnotationsToPeaks.R")
source("/data/Sox8_binding_partner_analysis/Functions/expand_ranges.R")

# Load the ArchR projects
ss4_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR/', force = FALSE)
ss8_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Sox8_Save_ArchR/', force = FALSE)

# Have to repeat the motif scanning on these ArchR projects since subsetting the peakset
ss4_Sox8_ArchRProj <- addMotifAnnotations(ss4_Sox8_ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)
ss8_Sox8_ArchRProj <- addMotifAnnotations(ss8_Sox8_ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)

saveArchRProject(ss4_Sox8_ArchRProj, "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR/")
saveArchRProject(ss8_Sox8_ArchRProj, "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Sox8_Save_ArchR/")

# Load the GRanges objects corresponding to accessibile peaks in Placodal, NC, and accessible across both
ss4_Sox8_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_Placodal_peaks.rds")
ss4_Sox8_NC_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_peaks.rds")
ss4_Sox8_NC_and_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_and_Placodal_peaks.rds")

ss8_Sox8_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_Placodal_peaks.rds")
ss8_Sox8_NC_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_peaks.rds")
ss8_Sox8_NC_and_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_and_Placodal_peaks.rds")


#### PREPARE GRANGES OBJECTS FOR MOTIF SCANNING. SEARCHING EITHER SIDE OF THE MOTIF, NOT OVERLAPPING THE MOTIF ITSELF ####

# First step, extract Sox8 motif positions from the ArchR projects. Comparing outputs from all archR projects to check that the Sox8 motif annotation is consistent
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR/', force = FALSE)

ss4_Sox8_motif_positions <- getPositions(ss4ArchRProj, name = "Motif", annoName = "SOX8")
ss4_Sox8_motif_positions_sub <- getPositions(ss4_Sox8_ArchRProj, name = "Motif", annoName = "SOX8")
ss8_Sox8_motif_positions <- getPositions(ss8ArchRProj, name = "Motif", annoName = "SOX8")
ss8_Sox8_motif_positions_sub <- getPositions(ss8_Sox8_ArchRProj, name = "Motif", annoName = "SOX8")

# Only need to use one of these GRanges going forwards, as they are all the same
Sox8_motif_positions <- as(ss4_Sox8_motif_positions$SOX8, "GRanges")

# Subset Sox8_motif_positions by overlap with placodal/NC peaksets
# Using a custom function to find overlaps in IRanges and pull metadata over from peaksets 
ss4_Placodal_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_Sox8_Placodal_peaks)
ss4_NC_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_Sox8_NC_peaks)
ss4_NC_and_Placodal_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_Sox8_NC_and_Placodal_peaks)

ss8_Placodal_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_Sox8_Placodal_peaks)
ss8_NC_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_Sox8_NC_peaks)
ss8_NC_and_Placodal_Sox8_Positions<- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_Sox8_NC_and_Placodal_peaks)

# Use a custom function to expand the search area upstream and downstream of the Sox8 binding motif
# Specify the expansion amount either side in basepairs. Choose whether or not to include the Sox8 motif itself.
# If you exclude the Sox8 motif, this creates a GRanges with double the number of rows, as separate IRanges are needed for the upstream and downstream regions for each Sox8 motif
ss4_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss4_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss4_Placodal_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss4_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)

ss4_NC_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss4_NC_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss4_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss4_NC_and_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)

ss4_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss8_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss4_NC_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss8_NC_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss4_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss8_NC_and_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
