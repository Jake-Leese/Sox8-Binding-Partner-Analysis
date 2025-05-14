# Script to read in meme SEA enrichment sites, match these back to the original reproducible peak sets, and export as bed files for enhancer annotation

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

unique(ss4_pPPR_Peaks$peakType)

# Combine and remove duplicates using reduce()
ss4_NC_all_Peaks <- GenomicRanges::reduce(c(ss4_dNC_Peaks, ss4_NC_Peaks))
ss4_PPR_all_Peaks <- GenomicRanges::reduce(c(ss4_pPPR_Peaks, ss4_aPPR_Peaks))
ss8_NC_all_Peaks <- GenomicRanges::reduce(c(ss8_dNC_Peaks, ss8_NC_Peaks))
ss8_PPR_all_Peaks <- GenomicRanges::reduce(c(ss8_pPPR_Peaks, ss8_aPPR_Peaks))

# Remove "chr" from seqnames
seqlevels(ss4_NC_all_Peaks) <- gsub("chr", "", seqlevels(ss4_NC_all_Peaks))
seqlevels(ss4_PPR_all_Peaks) <- gsub("chr", "", seqlevels(ss4_PPR_all_Peaks))
seqlevels(ss8_NC_all_Peaks) <- gsub("chr", "", seqlevels(ss8_NC_all_Peaks))
seqlevels(ss8_PPR_all_Peaks) <- gsub("chr", "", seqlevels(ss8_PPR_all_Peaks))

###########################################################################################
#################### Load SEA motif matches and map back to rep peaks #####################

# Loading in Sox8+ peak enrichment analysis results
PPR_File_path <- "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_Outputs/SEA_outputs/Reproducible_PPR_peaks/"

ss8_PPR_SEA_motif_sites_PPR_bg <- read.table(paste0(PPR_File_path, "ss8_PPR_rep_SOX8_peaks_enrichment/sites.tsv"))
ss8_PPR_SEA_motif_sites_shuffled_bg <- read.table(paste0(PPR_File_path, "shuffled_bg/ss8_PPR_rep_SOX8_peaks_enrichment/sites.tsv"))

ss8_PPR_SEA_DLX4_sites_PPR_bg <- ss8_PPR_SEA_motif_sites %>% filter(V1 == "Dlx4")
ss8_PPR_SEA_LMX1A_sites_PPR_bg <- ss8_PPR_SEA_motif_sites %>% filter(V1 == "LMX1A")
ss8_PPR_SEA_DLX4_sites_PPR_bg_gr <- makeGRangesFromDataFrame(ss8_PPR_SEA_DLX4_sites_PPR_bg,
                                                      seqnames.field = "V3",
                                                      start.field = "V4",
                                                      end.field = "V5",
                                                      strand.field = "V6",
                                                      keep.extra.columns = TRUE)  # Keep other columns if present
ss8_PPR_SEA_LMX1A_sites_PPR_bg_gr <- makeGRangesFromDataFrame(ss8_PPR_SEA_LMX1A_sites_PPR_bg,
                                                             seqnames.field = "V3",
                                                             start.field = "V4",
                                                             end.field = "V5",
                                                             strand.field = "V6",
                                                             keep.extra.columns = TRUE)  # Keep other columns if present

seqlevels(ss8_PPR_SEA_DLX4_sites_PPR_bg_gr) <- gsub("chr", "", seqlevels(ss8_PPR_SEA_DLX4_sites_PPR_bg_gr))
seqlevels(ss8_PPR_SEA_LMX1A_sites_PPR_bg_gr) <- gsub("chr", "", seqlevels(ss8_PPR_SEA_LMX1A_sites_PPR_bg_gr))

# Function to find overlaps and subset subsect by matches
Keep_overlaps <- function(subject, query){
  # Find overlaps
  overlapping_indices <- queryHits(findOverlaps(subject, query))
  # Subset the GRanges object X to keep only overlapping ranges
  subject_overlapping <- subject[overlapping_indices]
  # Print the result
  return(subject_overlapping)
}

# Returning GRanges of peaks that have SOX8 and enriched secondary motif
ss8_PPR_SOX8_DLX4_PPR_bg_peaks <- Keep_overlaps(ss8_PPR_all_Peaks, ss8_PPR_SEA_DLX4_sites_gr)
ss8_PPR_SOX8_LMX1A_PPR_bg_peaks <- Keep_overlaps(ss8_PPR_all_Peaks, ss8_PPR_SEA_LMX1A_sites_gr)

NC_bed_path_PPR_bg <- "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/NC/Sea_enriched/PPR_bg/"
PPR_bed_path_PPR_bg <- "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/PPR/Sea_enriched/PPR_bg/"
NC_bed_path_shuffled_bg <- "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/NC/Sea_enriched/shuffled_bg/"
PPR_bed_path_shuffled_bg <- "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/PPR/Sea_enriched/shuffled_bg/"





