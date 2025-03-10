# Script to scan for SOX8 binding motifs and subset our peaksets accordingly. 
# Outputs subsetted peaks centered around Sox8 motifs for SpaMo analysis

.libPaths("/R/libs/AT_ArcR_macs2")
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
library(JASPAR2020)
library(motifmatchr)
library(universalmotif)

source("/data/Sox8_binding_partner_analysis/src/Functions/MatchAnnotationsToPeaks.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/expand_ranges.R")

path <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/"

# Read in All NC and PPR marker peaks from ss4 and ss8, and combine them 

# Read in full peaksets
ss4_PPR_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss4_PPR_all_Peaks.gr.rds")
ss4_NC_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss4_NC_all_Peaks.gr.rds")
ss8_PPR_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss8_PPR_all_Peaks.gr.rds")
ss8_NC_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss8_NC_all_Peaks.gr.rds")

# Combine the GRanges objects and remove duplicates
NC_PPR_combined_gr <- unique(c(ss4_PPR_all_Peaks, ss4_NC_all_Peaks, ss8_PPR_all_Peaks, ss8_NC_all_Peaks))


#########################################################################################################
######################################## MOTIF SCANNING #################################################

# Define background nucleotide frequencies (based on full ArchR peakset)
ACTG_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Get the SOX8 motif from JASPAR2020
SOX8_tfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM", name = "SOX8"))[[1]]

# Scan for motif occurrences in the peakset
motif_hits <- matchMotifs(SOX8_tfm, NC_PPR_combined_gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                          out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
saveRDS(motif_hits, paste0(path, "Sox8_positions_all_peaks.gr.rds"))

# Subset original GRanges object to only those with a SOX8 binding motif
NC_PPR_combined_SOX8_gr_orig <- subsetByOverlaps(NC_PPR_combined_gr, motif_hits)


########################################################################################################
############################ Subsetting SOX8 peaks to NC, PPR, and shared ##############################

# Reading in full peaksets (not just SOX8 pos) 
ss4_PPR_only <- readRDS(paste0(path, "ss4_PPR_only.gr.rds"))
ss4_NC_only <- readRDS(paste0(path, "ss4_NC_only.gr.rds"))
ss4_PPR_NC_Shared <- readRDS(paste0(path, "ss4_PPR_and_NC_shared.gr.rds"))
ss8_PPR_only <- readRDS(paste0(path, "ss8_PPR_only.gr.rds"))
ss8_NC_only <- readRDS(paste0(path, "ss8_NC_only.gr.rds"))
ss8_PPR_NC_Shared <- readRDS(paste0(path, "ss8_PPR_and_NC_shared.gr.rds"))

# Subsetting GRanges by ranges that have a Sox8 motif
ss4_PPR_Sox8 <- ss4_PPR_only[ss4_PPR_only %over% NC_PPR_combined_SOX8_gr] #296 peaks
ss4_NC_Sox8 <- ss4_NC_only[ss4_NC_only %over% NC_PPR_combined_SOX8_gr] # 674 peaks
ss4_NC_and_PPR_Sox8 <- ss4_PPR_NC_Shared[ss4_PPR_NC_Shared %over% NC_PPR_combined_SOX8_gr] # 312 peaks
ss8_PPR_Sox8 <- ss8_PPR_only[ss8_PPR_only %over% NC_PPR_combined_SOX8_gr] # 1236 peaks
ss8_NC_Sox8 <- ss8_NC_only[ss8_NC_only %over% NC_PPR_combined_SOX8_gr] # 1593 peaks
ss8_NC_and_PPR_Sox8 <- ss8_PPR_NC_Shared[ss8_PPR_NC_Shared %over% NC_PPR_combined_SOX8_gr] # 438 peaks

# Saving Sox8 peaks
saveRDS(ss4_PPR_Sox8, paste0(path, "ss4_PPR_Sox8.gr.rds"))
saveRDS(ss4_NC_Sox8, paste0(path, "ss4_NC_Sox8.gr.rds"))
saveRDS(ss4_NC_and_PPR_Sox8, paste0(path, "ss4_NC_and_PPR_Sox8.gr.rds"))
saveRDS(ss8_PPR_Sox8, paste0(path, "ss8_PPR_Sox8.gr.rds"))
saveRDS(ss8_NC_Sox8, paste0(path, "ss8_NC_Sox8.gr.rds"))
saveRDS(ss8_NC_and_PPR_Sox8, paste0(path, "ss8_NC_and_PPR_Sox8.gr.rds"))

########################################################################################################
############################### Outputting peaksets for MEME analysis  #################################

#### PREPARE GRANGES CENTERED AROUND SOX8 MOTIF FOR DIFFERENT STAGE/CELL TYPES ####

# Load the GRanges objects corresponding to accessibile peaks in Placodal, NC, and accessible across both
ss4_PPR_Sox8 <- readRDS(paste0(path, "ss8_NC_and_PPR_Sox8.gr.rds"))        # 438 peaks
ss4_NC_Sox8 <- readRDS(paste0(path, "ss4_NC_Sox8.gr.rds"))                 # 675 peaks
ss4_NC_and_PPR_Sox8 <- readRDS(paste0(path, "ss4_NC_and_PPR_Sox8.gr.rds")) # 312 peaks

ss8_PPR_Sox8 <- readRDS(paste0(path, "ss8_PPR_Sox8.gr.rds"))               # 1236 peaks
ss8_NC_Sox8 <- readRDS(paste0(path, "ss8_NC_Sox8.gr.rds"))                 # 1593 peaks
ss8_NC_and_PPR_Sox8 <- readRDS(paste0(path, "ss8_NC_and_PPR_Sox8.gr.rds")) # 438 peaks

# Load the Sox8 motif positions
Sox8_motif_positions <- readRDS(paste0(path, "Sox8_positions_all_peaks.gr.rds"))

# Subset Sox8_motif_positions by overlap with placodal/NC peaksets
# Using a custom function to find overlaps in IRanges and pull metadata over from peaksets 
ss4_PPR_Sox8 <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_PPR_Sox8)
ss4_NC_Sox8 <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_NC_Sox8)
ss4_NC_and_PPR_Sox8 <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_NC_and_PPR_Sox8)

ss8_PPR_Sox8 <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_PPR_Sox8)
ss8_NC_Sox8 <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_NC_Sox8)
ss8_NC_and_PPR_Sox8 <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_NC_and_PPR_Sox8)


#
