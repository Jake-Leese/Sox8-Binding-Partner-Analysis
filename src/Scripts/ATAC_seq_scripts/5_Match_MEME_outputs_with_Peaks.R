# Script to import genome co-ordinates from MEME analysis and filter out which co-ordinates relate to which ArchR_Proj peaks along with the corresponding metadata
# This script assumes that the Sox8 binding motif is at the center of the co-ordinates

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
source("/data/Sox8_binding_partner_analysis/Functions/Match_GRanges_to_Metadata.R")

#### PREPARE GRANGES OF SOX8 MOTIF POSITIONS AND PULL OVER THEIR METADATA ####

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


#### IMPORT MEME ANALYSIS OUTPUTS FOR PEAKS OF INTEREST AND MATCH THEM BACK TO THE ORIGINAL PEAK METADATA ####

# Import the co-ordinates saved from SpaMo. In this example, for 3 different downstream (opposite strand) distances from the Sox8 motif
SOX4_2bp_Ds <- read.table("SpaMo_sequence_IDs/Round5_results/HH9_NC_TFs_NC_Peaks/SOX4/Ds_2.txt")
SOX4_3bp_Ds <- read.table("SpaMo_sequence_IDs/Round5_results/HH9_NC_TFs_NC_Peaks/SOX4/Ds_3.txt")
SOX4_4bp_Ds <- read.table("SpaMo_sequence_IDs/Round5_results/HH9_NC_TFs_NC_Peaks/SOX4/Ds_4.txt")

# Combine the co-ordinates and make into a GRanges object
Sox4_Ds_peaks <- c(SOX4_2bp_Ds[[1]], SOX4_3bp_Ds[[1]], SOX4_4bp_Ds[[1]])#

# Split the character vector into components (seqname, start, end)
Sox4_Ds_peaks_split <- strsplit(Sox4_Ds_peaks, "[:-]")

# Convert to a GRanges object
Sox4_Ds_peaks_gr <- GRanges(
  seqnames = sapply(Sox4_Ds_peaks_split, function(x) x[1]),             # Sequence names (chromosomes)
  ranges = IRanges(
    start = as.numeric(sapply(Sox4_Ds_peaks_split, function(x) x[2])),  # Start positions
    end = as.numeric(sapply(Sox4_Ds_peaks_split, function(x) x[3]))     # End positions
  )
)

Sox4_Ds_peaks_gr <- sort(Sox4_Ds_peaks_gr) # Sorts ranges in order

# This custom function finds overlaps with our output co-ordinates against the Sox8 motif positions and attaches the metadata from the 
# positions GRanges object, which is essentially the GRanges object for the corresponding ArchR peaks
Sox4_Ds_peaks_gr <- MatchAnnotationsToPeaks(Sox4_Ds_peaks_gr, ss8_NC_Sox8_Positions)

# Re-order GRanges for seqnames to be in order
mcols(Sox4_Ds_peaks_gr)$chr_number <- as.numeric(sub("chr", "", Sox4_Ds_peaks_gr$seqnames))
Sox4_Ds_peaks_gr <- Sox4_Ds_peaks_gr[order(mcols(Sox4_Ds_peaks_gr)$chr_number)] 



#### IMPORT PEAK2GENELINKS FOR OUR PEAKS OF INTEREST TO SEE WHAT GENES THEY MAY BE ASSOCIATED WITH ####

ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

# Retrieve the Peak2GeneLinks dataframe, which returns an ATAC peak (idxATAC) linked to a gene (idxRNA) based on co-ordinated accessibility and expression
ss8_p2g <- getPeak2GeneLinks(
  ss8ArchRProj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE)

# Get the ArchR peaks from the metadata stored in the above object. 
ArchR_peaks <- metadata(ss8_p2g)[[1]]

# Get the gene names from the metadata stored in the above object.
gene_names <- metadata(ss8_p2g)$geneSet$name

# Add the gene names to the p2g dataframe. sub(".*:", "", ) replaces the values of the vector up to the colon and replaces with "". It is just a method to remove the character string before the gene name
ss8_p2g$gene_names <- gene_names[ss8_p2g$idxRNA]

# Make a GRanges object from the idxATAC peaks and add the p2g dataframe as metadata columns

ArchR_peaks_subset <- ArchR_peaks[ss8_p2g$idxATAC]

mcols(ArchR_peaks_subset) <- ss8_p2g
ss8_p2g_gr <- sort(ArchR_peaks_subset)

saveRDS(ss8_p2g_gr, "Peaksets/ss8_Peaks2Genes.rds")

# Match the peaks2genes peaks with those of interest from out MEME analysis
ss8_p2g_gr <- readRDS("Peaksets/ss8_Peaks2Genes.rds")
Sox4_Ds_peaks_p2g <- Match_GRanges_to_Metadata(ss8_p2g_gr, Sox4_Ds_peaks_gr)  # 0 hits with peaks2Genes peak list

# Checking ss8 peak2genes against ss8_NC accessible peaks.
ss8_Sox8_NC_peaks_p2g <- sort(Match_GRanges_to_Metadata(ss8_p2g_gr, ss8_Sox8_NC_peaks))
ss8_Sox8_NC_peaks_p2g$gene_names    # 34 matches and only 31 linked genes, list doesn't contain commonly known NC genes, suggesting many could be missed


#### ALTERNATIVE: LOOK AT LIST OF NEAREST GENES (AS GIVEN IN METADATA) TO SOX4 ENRICHED PEAKS ####

ss8_Sox4_Sox8_nearest_genes <- Sox4_Ds_peaks_gr$nearestGene
Sox4_Ds_peaks_gr[Sox4_Ds_peaks_gr$nearestGene == "SNAI2"]

 
#### PLOT ARCHR BROWSER TRACKS FOR SOX4 PEAKS OF INTEREST THAT FALL IN THE RANGE OF THE P2G GRANGES ####

ss8_p2g_loops <- getPeak2GeneLinks(
  ss8ArchRProj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE)
ss8_p2g_loops <- ss8_p2g_loops$Peak2GeneLinks

# Add gene names to the GRanges
mcols(ss8_p2g_loops)$idxATAC <- ss8_p2g_gr_sorted$idxATAC
mcols(ss8_p2g_loops)$idxRNA <- ss8_p2g_gr_sorted$idxRNA
mcols(ss8_p2g_loops)$gene_names <- sub(".*:", "", gene_names[mcols(ss8_p2g_loops)$idxRNA])

Sox4_Ds_peaks_gr[subjectHits(findOverlaps(ss8_p2g_loops, Sox4_Ds_peaks_gr))]

MatchAnnotationsToPeaks(ss8_p2g_loops, Sox4_Ds_peaks_gr)$gene_names
Match_GRanges_to_Metadata(Sox4_Ds_peaks_gr, ss8_p2g_loops)
