# Script for peak calling peaks for broad_sc_helper_celltypes in ss4 and ss8 ArchR projects
# This is with the aim of identifying all open peaks in placodal and NC cell types, rather than just focussing on differential accessibility.
# The output should give us GRanges objects for all accessible Placodal and NC peaks for ss4 and ss8 stages. These can then be further subsetted into NC- or placode-only peaksets for downstream analysis
# Using alexthiery-schelper-archr_dev_macs2-schelper-0.3.5.img container

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
library(JASPAR2020)
library(BSgenome.Ggallus.UCSC.galGal6)
library(org.Gg.eg.db)
library(motifmatchr)
library(TFBSTools)


# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_celltype_peaks_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_celltype_peaks_Save-ArchR', force = FALSE)

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
  force = TRUE,
  genomeSize = 1230258557, # copied from Grace's paper, need to check this
)

ss8ArchRProj <- addReproduciblePeakSet(
  ArchRProj = ss8ArchRProj,
  groupBy = "transferred_scHelper_cell_type",
  pathToMacs2 = "/opt/conda/bin/macs2",
  force = TRUE,
  genomeSize = 1230258557, # copied from Grace's paper, need to check this
)


##############################################################################################
############# Read in reproducible peaksets for transferred_scHelper_cell_types ##############

# Read in cell type peaksets
ss4_pPPR_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/pPPR-reproduciblePeaks.gr.rds")
ss4_aPPR_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/aPPR-reproduciblePeaks.gr.rds")
ss4_dNC_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/dNC-reproduciblePeaks.gr.rds")
ss4_NC_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/NC-reproduciblePeaks.gr.rds")

ss8_pPPR_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/pPPR-reproduciblePeaks.gr.rds")
ss8_aPPR_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/aPPR-reproduciblePeaks.gr.rds")
ss8_dNC_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/dNC-reproduciblePeaks.gr.rds")
ss8_NC_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_celltype_peaks_Save-ArchR/PeakCalls/transferred_scHelper_cell_type/NC-reproduciblePeaks.gr.rds")

# Combine and remove duplicates using reduce()
ss4_NC_all_Peaks <- GenomicRanges::reduce(c(ss4_dNC_Peaks, ss4_NC_Peaks))
ss4_PPR_all_Peaks <- GenomicRanges::reduce(c(ss4_pPPR_Peaks, ss4_aPPR_Peaks))
ss8_NC_all_Peaks <- GenomicRanges::reduce(c(ss8_dNC_Peaks, ss8_NC_Peaks))
ss8_PPR_all_Peaks <- GenomicRanges::reduce(c(ss8_pPPR_Peaks, ss8_aPPR_Peaks))

# Confirming duplicates are limited to ranges, not seqnames and ranges
ss4_NC_all_Peaks[ranges(ss4_NC_all_Peaks) %in% ranges(ss4_NC_all_Peaks)[duplicated(ranges(ss4_NC_all_Peaks))]]
ss4_PPR_all_Peaks[ranges(ss4_PPR_all_Peaks) %in% ranges(ss4_PPR_all_Peaks)[duplicated(ranges(ss4_PPR_all_Peaks))]]
ss8_NC_all_Peaks[ranges(ss8_NC_all_Peaks) %in% ranges(ss8_NC_all_Peaks)[duplicated(ranges(ss8_NC_all_Peaks))]]
ss8_PPR_all_Peaks[ranges(ss8_PPR_all_Peaks) %in% ranges(ss8_PPR_all_Peaks)[duplicated(ranges(ss8_PPR_all_Peaks))]]


#########################################################################################################
######################################## MOTIF SCANNING #################################################

# Define background nucleotide frequencies (based on full ArchR peakset)
ACTG_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Get the SOX8 motif from JASPAR2020
SOX8_tfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM", name = "SOX8"))[[1]]

# Scan for motif occurrences in the peakset
ss4_NC_motif_hits <- matchMotifs(SOX8_tfm, ss4_NC_all_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                               out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss4_PPR_motif_hits <- matchMotifs(SOX8_tfm, ss4_PPR_all_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                 out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss4_pPPR_motif_hits <- matchMotifs(SOX8_tfm, ss4_pPPR_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                   out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss8_NC_motif_hits <- matchMotifs(SOX8_tfm, ss8_NC_all_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                 out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss8_PPR_motif_hits <- matchMotifs(SOX8_tfm, ss8_PPR_all_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                 out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss8_pPPR_motif_hits <- matchMotifs(SOX8_tfm, ss8_pPPR_Peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                   out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)

source("/data/Sox8_binding_partner_analysis/src/Functions/expand_ranges.R")

# Expand SOX8 motif hits +/- 100 bp either side
ss4_NC_motif_hits_100bp_surr <- expand_ranges(ss4_NC_motif_hits, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_PPR_motif_hits_100bp_surr <- expand_ranges(ss4_PPR_motif_hits, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_pPPR_motif_hits_100bp_surr <- expand_ranges(ss4_pPPR_motif_hits, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_motif_hits_100bp_surr <- expand_ranges(ss8_NC_motif_hits, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_PPR_motif_hits_100bp_surr <- expand_ranges(ss8_PPR_motif_hits, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_pPPR_motif_hits_100bp_surr <- expand_ranges(ss8_pPPR_motif_hits, expansion_amount = 100, KeepExistingRange = TRUE)


#### CREATE AND SAVE FASTA FILES FOR PEAKSETS AS INPUTS FOR MEME ####

GalGal6 <- BSgenome.Ggallus.UCSC.galGal6

# Get raw sequence using the GRanges objects
ss4_NC_motif_hits_100bp_seq <- getSeq(GalGal6, ss4_NC_motif_hits_100bp_surr)
ss4_PPR_motif_hits_100bp_seq <- getSeq(GalGal6, ss4_PPR_motif_hits_100bp_surr)
ss4_pPPR_motif_hits_100bp_seq <- getSeq(GalGal6, ss4_pPPR_motif_hits_100bp_surr)
ss8_NC_motif_hits_100bp_seq <- getSeq(GalGal6, ss8_NC_motif_hits_100bp_surr)
ss8_PPR_motif_hits_100bp_seq <- getSeq(GalGal6, ss8_PPR_motif_hits_100bp_surr)
ss8_pPPR_motif_hits_100bp_seq <- getSeq(GalGal6, ss8_pPPR_motif_hits_100bp_surr)

# Adding seqnames and ranges as unique sequence identifiers
names(ss4_NC_motif_hits_100bp_seq) <- paste0(ss4_NC_motif_hits_100bp_surr@seqnames, ":", ss4_NC_motif_hits_100bp_surr@ranges)
names(ss4_PPR_motif_hits_100bp_seq) <- paste0(ss4_PPR_motif_hits_100bp_surr@seqnames, ":", ss4_PPR_motif_hits_100bp_surr@ranges)
names(ss4_pPPR_motif_hits_100bp_seq) <- paste0(ss4_pPPR_motif_hits_100bp_surr@seqnames, ":", ss4_pPPR_motif_hits_100bp_surr@ranges)
names(ss8_NC_motif_hits_100bp_seq) <- paste0(ss8_NC_motif_hits_100bp_surr@seqnames, ":", ss8_NC_motif_hits_100bp_surr@ranges)
names(ss8_PPR_motif_hits_100bp_seq) <- paste0(ss8_PPR_motif_hits_100bp_surr@seqnames, ":", ss8_PPR_motif_hits_100bp_surr@ranges)
names(ss8_pPPR_motif_hits_100bp_seq) <- paste0(ss8_pPPR_motif_hits_100bp_surr@seqnames, ":", ss8_pPPR_motif_hits_100bp_surr@ranges)

# Write out sequences as .fasta files
path_MEME <- "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/fasta_files/"
writeXStringSet(ss4_NC_motif_hits_100bp_seq, file=paste0(path_MEME, "ss4_NC_rep_peaks_SOX8_100bp.fasta"))
writeXStringSet(ss4_PPR_motif_hits_100bp_seq, file=paste0(path_MEME, "ss4_PPR_rep_peaks_SOX8_100bp.fasta"))
writeXStringSet(ss4_pPPR_motif_hits_100bp_seq, file=paste0(path_MEME, "ss4_pPPR_rep_peaks_SOX8_100bp.fasta"))
writeXStringSet(ss8_NC_motif_hits_100bp_seq, file=paste0(path_MEME, "ss8_NC_rep_peaks_SOX8_100bp.fasta"))
writeXStringSet(ss8_PPR_motif_hits_100bp_seq, file=paste0(path_MEME, "ss8_PPR_rep_peaks_SOX8_100bp.fasta"))
writeXStringSet(ss8_pPPR_motif_hits_100bp_seq, file=paste0(path_MEME, "ss8_pPPR_rep_peaks_SOX8_100bp.fasta"))

