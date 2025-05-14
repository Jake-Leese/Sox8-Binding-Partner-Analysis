# Script for preparing PPR and NC differentially accessible peaks for Six1 SpaMo analysis

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
path_MEME <- "/data/Meme_suite/Six1_NC_v_Placodes/fasta_files/"

# Output directory
output_dir <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/Pairwise_peak_comparisons/"

# Read back in GRanges Lists
ss4_PPR_accessible_peaks <- readRDS(paste0(output_dir, "ss4/PPR_pairwise_peaks.gr.list.rds"))
ss4_NC_accessible_peaks <- readRDS(paste0(output_dir, "ss4/NC_pairwise_peaks.gr.list.rds"))
ss8_PPR_accessible_peaks <- readRDS(paste0(output_dir, "ss8/PPR_pairwise_peaks.gr.list.rds"))
ss8_NC_accessible_peaks <- readRDS(paste0(output_dir, "ss8/NC_pairwise_peaks.gr.list.rds"))

# Extract Differentially accessible peaks between NC and PPR for each stage
ss4_PPR_DE_peaks <- ss4_PPR_accessible_peaks$Placodal_vs_NC
ss4_NC_DE_peaks <- ss4_NC_accessible_peaks$NC_vs_Placodal
ss8_PPR_DE_peaks <- ss8_PPR_accessible_peaks$Placodal_vs_NC
ss8_NC_DE_peaks <- ss8_NC_accessible_peaks$NC_vs_Placodal


#########################################################################################################
######################################## MOTIF SCANNING #################################################


# Define background nucleotide frequencies (based on full ArchR peakset)
ACTG_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Get the SOX8 motif from JASPAR2020
SIX1_tfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM", name = "SIX1"))[[1]]

# Scan for motif occurrences in the peakset with different p. cutoffs
ss4_PPR_SIX1_motif_hits_e05 <- matchMotifs(SIX1_tfm, ss4_PPR_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                           out = "positions", p.cutoff = 1e-05, bg = ACTG_freqs)
ss4_PPR_SIX1_motif_hits_e04 <- matchMotifs(SIX1_tfm, ss4_PPR_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                               out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss4_NC_SIX1_motif_hits_e05 <- matchMotifs(SIX1_tfm, ss4_NC_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                           out = "positions", p.cutoff = 1e-05, bg = ACTG_freqs)
ss4_NC_SIX1_motif_hits_e04 <- matchMotifs(SIX1_tfm, ss4_NC_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                           out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss8_PPR_SIX1_motif_hits_e05 <- matchMotifs(SIX1_tfm, ss8_PPR_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                           out = "positions", p.cutoff = 1e-05, bg = ACTG_freqs)
ss8_PPR_SIX1_motif_hits_e04 <- matchMotifs(SIX1_tfm, ss8_PPR_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                           out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)
ss8_NC_SIX1_motif_hits_e05 <- matchMotifs(SIX1_tfm, ss8_NC_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                          out = "positions", p.cutoff = 1e-05, bg = ACTG_freqs)
ss8_NC_SIX1_motif_hits_e04 <- matchMotifs(SIX1_tfm, ss8_NC_DE_peaks, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                          out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)

# Expand GRanges (+/- 100 bp) around SIX1 motif positions
ss4_PPR_SIX1_motif_hits_e05_100bp_surr <- expand_ranges(ss4_PPR_SIX1_motif_hits_e05, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_PPR_SIX1_motif_hits_e04_100bp_surr <- expand_ranges(ss4_PPR_SIX1_motif_hits_e04, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_SIX1_motif_hits_e05_100bp_surr <- expand_ranges(ss4_NC_SIX1_motif_hits_e05, expansion_amount = 100, KeepExistingRange = TRUE)
ss4_NC_SIX1_motif_hits_e04_100bp_surr <- expand_ranges(ss4_NC_SIX1_motif_hits_e04, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_PPR_SIX1_motif_hits_e05_100bp_surr <- expand_ranges(ss8_PPR_SIX1_motif_hits_e05, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_PPR_SIX1_motif_hits_e04_100bp_surr <- expand_ranges(ss8_PPR_SIX1_motif_hits_e04, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_SIX1_motif_hits_e05_100bp_surr <- expand_ranges(ss8_NC_SIX1_motif_hits_e05, expansion_amount = 100, KeepExistingRange = TRUE)
ss8_NC_SIX1_motif_hits_e04_100bp_surr <- expand_ranges(ss8_NC_SIX1_motif_hits_e04, expansion_amount = 100, KeepExistingRange = TRUE)

# Convert and save as fasta files
prepare_fasta <- function(Genome, GRanges){
  fasta <- getSeq(Genome, GRanges)
  names(fasta) <- paste0(GRanges@seqnames, ":", GRanges@ranges)
  return(fasta)
}

GalGal6 <- BSgenome.Ggallus.UCSC.galGal6

ss4_PPR_SIX1_motif_hits_e05_100bp_seq <- prepare_fasta(GalGal6, ss4_PPR_SIX1_motif_hits_e05_100bp_surr)
ss4_PPR_SIX1_motif_hits_e04_100bp_seq <- prepare_fasta(GalGal6, ss4_PPR_SIX1_motif_hits_e04_100bp_surr)
ss4_NC_SIX1_motif_hits_e05_100bp_seq <- prepare_fasta(GalGal6, ss4_NC_SIX1_motif_hits_e05_100bp_surr)
ss4_NC_SIX1_motif_hits_e04_100bp_seq <- prepare_fasta(GalGal6, ss4_NC_SIX1_motif_hits_e04_100bp_surr)
ss8_PPR_SIX1_motif_hits_e05_100bp_seq <- prepare_fasta(GalGal6, ss8_PPR_SIX1_motif_hits_e05_100bp_surr)
ss8_PPR_SIX1_motif_hits_e04_100bp_seq <- prepare_fasta(GalGal6, ss8_PPR_SIX1_motif_hits_e04_100bp_surr)
ss8_NC_SIX1_motif_hits_e05_100bp_seq <- prepare_fasta(GalGal6, ss8_NC_SIX1_motif_hits_e05_100bp_surr)
ss8_NC_SIX1_motif_hits_e04_100bp_seq <- prepare_fasta(GalGal6, ss8_NC_SIX1_motif_hits_e04_100bp_surr)

# Write out sequences as .fasta files
writeXStringSet(ss4_PPR_SIX1_motif_hits_e05_100bp_seq, file=paste0(path_MEME, "ss4_PPR_SIX1_e05_100bp_surr.fasta"))
writeXStringSet(ss4_PPR_SIX1_motif_hits_e04_100bp_seq, file=paste0(path_MEME, "ss4_PPR_SIX1_e04_100bp_surr.fasta"))
writeXStringSet(ss4_NC_SIX1_motif_hits_e05_100bp_seq, file=paste0(path_MEME, "ss4_NC_SIX1_e05_100bp_surr.fasta"))
writeXStringSet(ss4_NC_SIX1_motif_hits_e04_100bp_seq, file=paste0(path_MEME, "ss4_NC_SIX1_e04_100bp_surr.fasta"))
writeXStringSet(ss8_PPR_SIX1_motif_hits_e05_100bp_seq, file=paste0(path_MEME, "ss8_PPR_SIX1_e05_100bp_surr.fasta"))
writeXStringSet(ss8_PPR_SIX1_motif_hits_e04_100bp_seq, file=paste0(path_MEME, "ss8_PPR_SIX1_e04_100bp_surr.fasta"))
writeXStringSet(ss8_NC_SIX1_motif_hits_e05_100bp_seq, file=paste0(path_MEME, "ss8_NC_SIX1_e05_100bp_surr.fasta"))
writeXStringSet(ss8_NC_SIX1_motif_hits_e04_100bp_seq, file=paste0(path_MEME, "ss8_NC_SIX1_e04_100bp_surr.fasta"))


# Preparing meme motif_matrix files
# Prepare motif lists
motifList_vert <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# rename each motif from their ID to their TF name
name_vector_vert <- c()
for (i in 1:length(motifList_vert)){
  name <- name(motifList_vert[[i]])
  name_vector_vert <- c(name_vector_vert, name)
}
names(motifList_vert) <- name_vector_vert

# Extract SIX1 and PAX7 motif matrices to use as primary and secondary motifs, respectively, in SpaMo analysis
SIX1_motif <- motifList_vert$SIX1
PAX7_motif <- motifList_vert$PAX7

# Full scATAC peakset background ACTG_freqs = 0.2487819, 0.2512324, 0.2512223, 0.2487634
ACGT_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Write JASPAR Motif matrices to minimal meme format
write_meme(SIX1_motif, "/data/Meme_suite/Six1_NC_v_Placodes/motifs/Six1.txt", version = 5, ACGT_freqs, overwrite = TRUE)

write_meme(PAX7_motif, "/data/Meme_suite/Six1_NC_v_Placodes/motifs/Pax7.txt", version = 5, ACGT_freqs, overwrite = TRUE)


