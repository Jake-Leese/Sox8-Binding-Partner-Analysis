# Script to write out fasta files for MEME suite motif enrichment/scanning #

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
library(JASPAR2024)
library(universalmotif)

source("/data/Sox8_binding_partner_analysis/Functions/MatchAnnotationsToPeaks.R")
source("/data/Sox8_binding_partner_analysis/Functions/expand_ranges.R")

#### PREPARE GRANGES CENTERED AROUND SOX8 MOTIF FOR DIFFERENT STAGE/CELL TYPES ####

# Load the GRanges objects corresponding to accessibile peaks in Placodal, NC, and accessible across both
ss4_Sox8_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_Placodal_peaks.rds")
ss4_Sox8_NC_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_peaks.rds")
ss4_Sox8_NC_and_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_and_Placodal_peaks.rds")

ss8_Sox8_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_Placodal_peaks.rds")
ss8_Sox8_NC_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_peaks.rds")
ss8_Sox8_NC_and_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_and_Placodal_peaks.rds")

# First step, extract all sox8 binding motif co-ordinates from an ArchR project
# ss4_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR/', force = FALSE)
# Sox8_motif_positions <- getPositions(ss4_Sox8_ArchRProj, name = "Motif", annoName = "SOX8")
# Sox8_motif_positions <- as(Sox8_motif_positions$SOX8, "GRanges")
# saveRDS(Sox8_motif_positions, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Sox8_motif_positions.rds")
Sox8_motif_positions <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Sox8_motif_positions.rds")

# Subset Sox8_motif_positions by overlap with placodal/NC peaksets
# Using a custom function to find overlaps in IRanges and pull metadata over from peaksets 
ss4_Placodal_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_Sox8_Placodal_peaks)
ss4_NC_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_Sox8_NC_peaks)
ss4_NC_and_Placodal_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss4_Sox8_NC_and_Placodal_peaks)

ss8_Placodal_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_Sox8_Placodal_peaks)
ss8_NC_Sox8_Positions <- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_Sox8_NC_peaks)
ss8_NC_and_Placodal_Sox8_Positions<- MatchAnnotationsToPeaks(Sox8_motif_positions, ss8_Sox8_NC_and_Placodal_peaks)

# Expand GRanges (+/- 250 bp) around Sox8 motif positions
ss4_Placodal_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss4_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)
ss4_NC_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss4_NC_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)
ss4_NC_and_Placodal_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss4_NC_and_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)

ss8_Placodal_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss8_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)
ss8_NC_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss8_NC_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)
ss8_NC_and_Placodal_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss8_NC_and_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)

# Expand GRanges (+/- 100 bp) around Sox8 motif positions
ss4_Placodal_Sox8_100bp_surrounding_incl_motif <- expand_ranges(ss4_Placodal_Sox8_Positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_Sox8_100bp_surrounding_incl_motif <- expand_ranges(ss4_NC_Sox8_Positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif <- expand_ranges(ss4_NC_and_Placodal_Sox8_Positions, expansion_amount = 100, KeepExistingRange = TRUE)

ss8_Placodal_Sox8_100bp_surrounding_incl_motif <- expand_ranges(ss8_Placodal_Sox8_Positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_Sox8_100bp_surrounding_incl_motif <- expand_ranges(ss8_NC_Sox8_Positions, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif <- expand_ranges(ss8_NC_and_Placodal_Sox8_Positions, expansion_amount = 100, KeepExistingRange = TRUE)

#### CREATE AND SAVE FASTA FILES AS INPUTS FOR MEME ####

GalGal6 <- BSgenome.Ggallus.UCSC.galGal6

# Get raw sequence using the GRanges objects
ss4_Placodal_Sox8_seq <- getSeq(GalGal6, ss4_Placodal_Sox8_100bp_surrounding_incl_motif)
ss4_NC_Sox8_seq <- getSeq(GalGal6, ss4_NC_Sox8_100bp_surrounding_incl_motif)
ss4_NC_and_Placodal_Sox8_seq <- getSeq(GalGal6, ss4_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif)
ss8_Placodal_Sox8_seq <- getSeq(GalGal6, ss8_Placodal_Sox8_100bp_surrounding_incl_motif)
ss8_NC_Sox8_seq <- getSeq(GalGal6, ss8_NC_Sox8_100bp_surrounding_incl_motif)
ss8_NC_and_Placodal_Sox8_seq <- getSeq(GalGal6, ss8_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif)

# Adding seqnames and ranges as unique sequence identifiers
names(ss4_Placodal_Sox8_seq) <- paste0(ss4_Placodal_Sox8_100bp_surrounding_incl_motif@seqnames, ":", ss4_Placodal_Sox8_100bp_surrounding_incl_motif@ranges)
names(ss4_NC_Sox8_seq) <- paste0(ss4_NC_Sox8_100bp_surrounding_incl_motif@seqnames, ":", ss4_NC_Sox8_100bp_surrounding_incl_motif@ranges)
names(ss4_NC_and_Placodal_Sox8_seq) <- paste0(ss4_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif@seqnames, ":", ss4_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif@ranges)
names(ss8_Placodal_Sox8_seq) <- paste0(ss8_Placodal_Sox8_100bp_surrounding_incl_motif@seqnames, ":", ss8_Placodal_Sox8_100bp_surrounding_incl_motif@ranges)
names(ss8_NC_Sox8_seq) <- paste0(ss8_NC_Sox8_100bp_surrounding_incl_motif@seqnames, ":", ss8_NC_Sox8_100bp_surrounding_incl_motif@ranges)
names(ss8_NC_and_Placodal_Sox8_seq) <- paste0(ss8_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif@seqnames, ":", ss8_NC_and_Placodal_Sox8_100bp_surrounding_incl_motif@ranges)

# Write out sequences as .fasta files
writeXStringSet(ss4_Placodal_Sox8_seq, file="/data/Sox8_binding_partner_analysis/meme_suite/fasta_files/ss4_Placodal_Sox8.fasta")
writeXStringSet(ss4_NC_Sox8_seq, file="/data/Sox8_binding_partner_analysis/meme_suite/fasta_files/ss4_NC_Sox8.fasta")
writeXStringSet(ss4_NC_and_Placodal_Sox8_seq, file="/data/Sox8_binding_partner_analysis/meme_suite/fasta_files/ss4_NC_and_Placodal_Sox8.fasta")
writeXStringSet(ss8_Placodal_Sox8_seq, file="/data/Sox8_binding_partner_analysis/meme_suite/fasta_files/ss8_Placodal_Sox8.fasta")
writeXStringSet(ss8_NC_Sox8_seq, file="/data/Sox8_binding_partner_analysis/meme_suite/fasta_files/ss8_NC_Sox8.fasta")
writeXStringSet(ss8_NC_and_Placodal_Sox8_seq, file="/data/Sox8_binding_partner_analysis/meme_suite/fasta_files/ss8_NC_and_Placodal_Sox8.fasta")

#### PREPARE MOTIF LISTS ####

Jaspar2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(Jaspar2024))
motifList_vert <- getMatrixSet(x = sq24, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PFM"))

# rename each motif from their ID to their TF name
name_vector_vert <- c()
for (i in 1:length(motifList_vert)){
  name <- name(motifList_vert[[i]])
  name_vector_vert <- c(name_vector_vert, name)
}
names(motifList_vert) <- name_vector_vert

# Read in SOX8 co-expression heatmap clusters and extract TFs based on their co-expression patterns
HH7_HH8_Sox8_Coexpression_TFs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/Heatmap_clusters/Sox8_coexpression_heatmap_clusters_HH7_HH8_min_0.2_20C.csv")
HH7_HH8_Sox8_Coexpression_TFs$X <- NULL

HH7_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH7 == "Placodal",]$Gene
HH7_NC_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH7 == "NC",]$Gene
HH7_NC_and_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH7 == "Both",]$Gene
HH8_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH8 == "Placodal",]$Gene
HH8_NC_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH8 == "NC",]$Gene
HH8_NC_and_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH8 == "Both",]$Gene

# subset motif list based on motif lists of interest
HH7_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_Placodal_TFs)]
HH7_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_NC_TFs)]
HH7_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_NC_and_Placodal_TFs)]
HH8_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_Placodal_TFs)]
HH8_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_NC_TFs)]
HH8_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_NC_and_Placodal_TFs)]

# Extract SOX8 motif matrix to use as primary motif in SpaMo analysis
Sox8_motif <- motifList_vert$SOX8

# Full scATAC peakset background ACTG_freqs = 0.2487819, 0.2512324, 0.2512223, 0.2487634
ACGT_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Write JASPAR Motif matrices to minimal meme format

write_meme(Sox8_motif, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/Sox8.txt", version = 5, ACGT_freqs)

write_meme(HH7_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH7_Placodal_motifs.txt", version = 5, ACGT_freqs)
write_meme(HH7_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH7_NC_motifs.txt", version = 5, ACGT_freqs)
write_meme(HH7_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH7_NC_and_Placodal_motifs.txt", version = 5, ACGT_freqs)
write_meme(HH8_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_Placodal_motifs.txt", version = 5, ACGT_freqs)
write_meme(HH8_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_NC_motifs.txt", version = 5, ACGT_freqs)
write_meme(HH8_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_NC_and_Placodal_motifs.txt", version = 5, ACGT_freqs)
