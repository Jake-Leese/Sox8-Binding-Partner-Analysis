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
library(motifmatchr)
#### 

source("/data/Sox8_binding_partner_analysis/Functions/MatchAnnotationsToPeaks.R")
source("/data/Sox8_binding_partner_analysis/Functions/expand_ranges.R")

# Load the ArchR projects
ss4_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR/', force = FALSE)
# ss8_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Sox8_Save_ArchR/', force = FALSE)

# Have to repeat the motif scanning on these ArchR projects since subsetting the peakset
# ss4_Sox8_ArchRProj <- addMotifAnnotations(ss4_Sox8_ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)
# ss8_Sox8_ArchRProj <- addMotifAnnotations(ss8_Sox8_ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)

# saveArchRProject(ss4_Sox8_ArchRProj, "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR/")
# saveArchRProject(ss8_Sox8_ArchRProj, "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Sox8_Save_ArchR/")

# Load the GRanges objects corresponding to accessibile peaks in Placodal, NC, and accessible across both
ss4_Sox8_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_Placodal_peaks.rds")
ss4_Sox8_NC_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_peaks.rds")
ss4_Sox8_NC_and_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_and_Placodal_peaks.rds")

ss8_Sox8_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_Placodal_peaks.rds")
ss8_Sox8_NC_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_peaks.rds")
ss8_Sox8_NC_and_Placodal_peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_and_Placodal_peaks.rds")


#### A. PREPARE GRANGES OBJECTS FOR MOTIF SCANNING. SEARCHING EITHER SIDE OF THE MOTIF, NOT OVERLAPPING THE MOTIF ITSELF ####

# First step, extract Sox8 motif positions from the ArchR projects. Comparing outputs from all archR projects to check that the Sox8 motif annotation is consistent
# ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR/', force = FALSE)
# ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR/', force = FALSE)

# ss4_Sox8_motif_positions_og <- getPositions(ss4ArchRProj, name = "Motif", annoName = "SOX8")
ss4_Sox8_motif_positions <- getPositions(ss4_Sox8_ArchRProj, name = "Motif", annoName = "SOX8")
# ss8_Sox8_motif_positions_og <- getPositions(ss8ArchRProj, name = "Motif", annoName = "SOX8")
# ss8_Sox8_motif_positions <- getPositions(ss8_Sox8_ArchRProj, name = "Motif", annoName = "SOX8")

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
# ss4_Placodal_Sox8_250bp_surrounding_incl_motif <- expand_ranges(ss4_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = TRUE)
ss4_NC_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss4_NC_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss4_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss4_NC_and_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)

ss8_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss8_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss8_NC_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss8_NC_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)
ss8_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif <- expand_ranges(ss8_NC_and_Placodal_Sox8_Positions, expansion_amount = 250, KeepExistingRange = FALSE)


#### B. PREPARE MOTIF MATRIX LISTS FOR SCANNING ####

# Make TF lists based off hierarchical clustering of HH7 and HH8 NC and Placodal Sox8 co-expression values from scRNA-seq
# The following dataframe corresponds to the heatmap "Sox8_coexpression_heatmap_HH7_HH8_0.2_min_20C.png" made specifying 20 clusters for high cluster resolution. See "RNA_seq_scripts/5_Clustering_HH7_and_HH8_only.R"
# The HH7 and HH8 columns reflect manual annotation of Sox8 coexpression patterns based off this heatmap

HH7_HH8_Sox8_Coexpression_TFs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/Heatmap_clusters/Sox8_coexpression_heatmap_clusters_HH7_HH8_min_0.2_20C.csv")
HH7_HH8_Sox8_Coexpression_TFs$X <- NULL

HH7_NC_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH7 == "NC",]$Gene
HH7_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH7 == "Placodal",]$Gene
HH7_NC_and_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH7 == "Both",]$Gene

HH8_NC_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH8 == "NC",]$Gene
HH8_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH8 == "Placodal",]$Gene
HH8_NC_and_Placodal_TFs <- HH7_HH8_Sox8_Coexpression_TFs[HH7_HH8_Sox8_Coexpression_TFs$HH8 == "Both",]$Gene

# Prepare the JASPAR2024 Pwm motif matrix and subset by TF lists above
Jaspar2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(Jaspar2024))
motifList_vert <- getMatrixSet(x = sq24, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# rename each motif from their ID to their TF name
name_vector_vert <- c()
for (i in 1:length(motifList_vert)){
  name <- name(motifList_vert[[i]])
  name_vector_vert <- c(name_vector_vert, name)
}
names(motifList_vert) <- name_vector_vert

# subset motif list based on motif lists of interest
HH7_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_NC_TFs)]
HH7_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_Placodal_TFs)]
HH7_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_NC_and_Placodal_TFs)]

HH8_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_NC_TFs)]
HH8_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_Placodal_TFs)]
HH8_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_NC_and_Placodal_TFs)]


#### C. SCAN SOX8 SURROUNDING SEQUENCES WITH TFs OF INTEREST ####

# Get peakset from ss8_ArchR project. 
Full_ss4_peakSet <- getPeakSet(ArchRProj = ss4ArchRProj)

# Get the ACGT frequencies from the full ss8 peakset
freqs <- alphabetFrequency(getSeq(BSgenome.Ggallus.UCSC.galGal6, Full_ss4_peakSet))
ACGT_freqs <- c(sum(freqs[,"A"])/sum(freqs),
                sum(freqs[,"C"])/sum(freqs),
                sum(freqs[,"G"])/sum(freqs),
                sum(freqs[,"T"])/sum(freqs))
# Original peakset bg ACTG_freqs = 0.2487819, 0.2512324, 0.2512223, 0.2487634

# Scan GRanges object for motifs
HH7_NC_Peaks_NC_motif_matches <- matchMotifs(
  HH7_NC_MotifList_vert, 
  ss4_NC_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH7_NC_Peaks_NC_and_Placodal_motif_matches <-  matchMotifs(
  HH7_NC_and_Placodal_MotifList_vert, 
  ss4_NC_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH7_NC_and_Placodal_Peaks_NC_motif_matches <- matchMotifs(
  HH7_NC_MotifList_vert, 
  ss4_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH7_NC_and_Placodal_Peaks_NC_and_Placodal_motif_matches <- matchMotifs(
  HH7_NC_and_Placodal_MotifList_vert, 
  ss4_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH7_NC_and_Placodal_Peaks_Placodal_motif_matches <- matchMotifs(
  HH7_Placodal_MotifList_vert, 
  ss4_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH7_Placodal_Peaks_NC_and_Placodal_motif_matches <- matchMotifs(
  HH7_NC_and_Placodal_MotifList_vert, 
  ss4_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH7_Placodal_Peaks_Placodal_motif_matches <- matchMotifs(
  HH7_Placodal_MotifList_vert, 
  ss4_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

colSums(assays(HH7_NC_Peaks_NC_motif_matches)$motifMatches)
colSums(assays(HH7_NC_Peaks_NC_and_Placodal_motif_matches)$motifMatches)
colSums(assays(HH7_NC_and_Placodal_Peaks_NC_motif_matches)$motifMatches)
colSums(assays(HH7_NC_and_Placodal_Peaks_NC_and_Placodal_motif_matches)$motifMatches)
colSums(assays(HH7_NC_and_Placodal_Peaks_Placodal_motif_matches)$motifMatches)
colSums(assays(HH7_Placodal_Peaks_NC_and_Placodal_motif_matches)$motifMatches)
colSums(assays(HH7_Placodal_Peaks_Placodal_motif_matches)$motifMatches)

# Repeat above for HH8
HH8_NC_Peaks_NC_motif_matches <- matchMotifs(
  HH8_NC_MotifList_vert, 
  ss8_NC_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH8_NC_Peaks_NC_and_Placodal_motif_matches <-  matchMotifs(
  HH8_NC_and_Placodal_MotifList_vert, 
  ss8_NC_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH8_NC_and_Placodal_Peaks_NC_motif_matches <- matchMotifs(
  HH8_NC_MotifList_vert, 
  ss8_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH8_NC_and_Placodal_Peaks_NC_and_Placodal_motif_matches <- matchMotifs(
  HH8_NC_and_Placodal_MotifList_vert, 
  ss8_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH8_NC_and_Placodal_Peaks_Placodal_motif_matches <- matchMotifs(
  HH8_Placodal_MotifList_vert, 
  ss8_NC_and_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH8_Placodal_Peaks_NC_and_Placodal_motif_matches <- matchMotifs(
  HH8_NC_and_Placodal_MotifList_vert, 
  ss8_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

HH8_Placodal_Peaks_Placodal_motif_matches <- matchMotifs(
  HH8_Placodal_MotifList_vert, 
  ss8_Placodal_Sox8_250bp_surrounding_excl_motif, 
  genome = "BSgenome.Ggallus.UCSC.galGal6",
  bg = ACGT_freqs,
  out = "scores", 
  p.cutoff = 1e-05)

colSums(assays(HH8_NC_Peaks_NC_motif_matches)$motifMatches)
colSums(assays(HH8_NC_Peaks_NC_and_Placodal_motif_matches)$motifMatches)
colSums(assays(HH8_NC_and_Placodal_Peaks_NC_motif_matches)$motifMatches)
colSums(assays(HH8_NC_and_Placodal_Peaks_NC_and_Placodal_motif_matches)$motifMatches)
colSums(assays(HH8_NC_and_Placodal_Peaks_Placodal_motif_matches)$motifMatches)
colSums(assays(HH8_Placodal_Peaks_NC_and_Placodal_motif_matches)$motifMatches)
colSums(assays(HH8_Placodal_Peaks_Placodal_motif_matches)$motifMatches)
