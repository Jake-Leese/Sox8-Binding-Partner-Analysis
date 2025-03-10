# Script for peak calling peaks for broad_sc_helper_celltypes in ss4 and ss8 ArchR projects
# This is with the aim of identifying all open peaks in placodal and NC cell types, rather than just focussing on differential accessibility.
# The output should give us GRanges objects for all accessible Placodal and NC peaks for ss4 and ss8 stages. These can then be further subsetted into NC- or placode-only peaksets for downstream analysis

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
library(parallel)
library(clustree)
library(ComplexHeatmap)
library(BSgenome.Ggallus.UCSC.galGal6)
library(scHelper)
library(ggrepel)

# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

##############################################################################################
############################## Generating pseudo-replicates ##################################

# Make pseudo replicates (only in the case that group coverages can't be found)
ss4_pseudo_replicates <- addGroupCoverages(ss4ArchRProj, groupBy = "transferred_scHelper_cell_type", returnGroups = TRUE, force = TRUE)
ss8_pseudo_replicates <- addGroupCoverages(ss8ArchRProj, groupBy = "transferred_scHelper_cell_type", returnGroups = TRUE, force = TRUE)

# Plot table to see which samples and groups the pseudo replicate cells come from 

png(paste0('/data/Sox8_binding_partner_analysis/scATACseq_objects/plots/pseudoreplicates/ss4_pseudoreplicate_cell_counts_per_sample_table.png'), height = 40, width = 30, units = 'cm', res = 400)
ArchR_PseudoreplicateCounts(ss4ArchRProj, ss4_pseudo_replicates, group_by = "Sample")
graphics.off()
png(paste0('/data/Sox8_binding_partner_analysis/scATACseq_objects/plots/pseudoreplicates/ss8_pseudoreplicate_cell_counts_per_sample_table.png'), height = 40, width = 30, units = 'cm', res = 400)
ArchR_PseudoreplicateCounts(ss8ArchRProj, ss8_pseudo_replicates, group_by = "Sample")
graphics.off()

png(paste0('/data/Sox8_binding_partner_analysis/scATACseq_objects/plots/pseudoreplicates/ss4pseudoreplicate_cell_counts_per_group_table.png'), height = 40, width = 70, units = 'cm', res = 400)
ArchR_PseudoreplicateCounts(ss4ArchRProj, ss4_pseudo_replicates, group_by = "transferred_scHelper_cell_type")
graphics.off()
png(paste0('/data/Sox8_binding_partner_analysis/scATACseq_objects/plots/pseudoreplicates/ss8pseudoreplicate_cell_counts_per_group_table.png'), height = 40, width = 70, units = 'cm', res = 400)
ArchR_PseudoreplicateCounts(ss8ArchRProj, ss8_pseudo_replicates, group_by = "transferred_scHelper_cell_type")
graphics.off()

#####  Make actual pseudo-replicates for peak calling:
ss4ArchRProj <- addGroupCoverages(ss4ArchRProj, groupBy = "transferred_scHelper_cell_type", returnGroups = FALSE, force = TRUE)
ss8ArchRProj <- addGroupCoverages(ss8ArchRProj, groupBy = "transferred_scHelper_cell_type", returnGroups = FALSE, force = TRUE)


##############################################################################################
############################## Call peaks on pseudo-replicates ###############################

ss4ArchRProj <- addReproduciblePeakSet(
  ArchRProj = ss4ArchRProj,
  groupBy = "transferred_scHelper_cell_type",
  pathToMacs2 = "/opt/conda/bin/macs2",
  reproducibility = "(n+1)/2", # This means that the majority of samples (pseudobulk replicates for each cell type) must have a peak call at the given locus
  force = TRUE,
  genomeSize = 1230258557, # copied from Grace's paper, need to check this
)

ss8ArchRProj <- addReproduciblePeakSet(
  ArchRProj = ss8ArchRProj,
  groupBy = "transferred_scHelper_cell_type",
  pathToMacs2 = "/opt/conda/bin/macs2",
  reproducibility = "(n+1)/2",
  force = TRUE,
  genomeSize = 1230258557, # copied from Grace's paper, need to check this
)


##############################################################################################
############# Read in reproducible peaksets for transferred_scHelper_cell_types ##############

# Read in cell type peaksets
ss4_pPPR_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/pPPR-reproduciblePeaks.gr.rds")
ss4_aPPR_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/aPPR-reproduciblePeaks.gr.rds")
ss4_dNC_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/dNC-reproduciblePeaks.gr.rds")
ss4_NC_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/NC-reproduciblePeaks.gr.rds")

# Combine and remove duplicates using reduce()
ss4_NC_all_Peaks <- GenomicRanges::reduce(c(ss4_dNC_Peaks, ss4_NC_Peaks))
all(ss4_dNC_Peaks %over% ss4_NC_all_Peaks) & all(ss4_NC_Peaks %over% ss4_NC_all_Peaks) # Should return TRUE if all original ranges are accounted for
ss4_NC_all_Peaks[which(duplicated(ranges(ss4_NC_all_Peaks)))] # Should read FALSE if there are no duplicates

anyDuplicated(ss4_PPR_all_Peaks)

ss4_pPPR_Peaks_ordered <- ss4_pPPR_Peaks[order(ss4_pPPR_Peaks)]

ss4_peak_matrix <- getMatrixFromProject(ss4ArchRProj, useMatrix = "PeakMatrix")
ss4_peak_counts <- assays(ss4_peak_matrix)$PeakMatrix

