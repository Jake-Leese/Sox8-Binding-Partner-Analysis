# Script to check overlap of peaks with motif of interest against given genome co-ordinates, in this case accessibility modules

.libPaths("/R/libs/AT_ArcR_macs2")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(ArchR)
library(tidyverse)
library(BSgenome.Ggallus.UCSC.galGal6)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(universalmotif)

source("/data/Sox8_binding_partner_analysis/src/Functions/convert_coordinates_to_GRanges.R")

##############################################################################################################################
#################################### Extracting peaks with our motif of interest #############################################

# loading ArchR project
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

# Manually read in motif matrix from ArchR project
motif_annotations <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR/Annotations/Motif-Matches-In-Peaks.rds")

# Extract the sparse matrix from the ranged experiment object
motif_annotations_matrix <- assays(motif_annotations)$matches
rownames(motif_annotations_matrix) <- rowData(motif_annotations)$name
colnames(motif_annotations_matrix) <- colnames(motif_annotations)
dim(motif_annotations_matrix)

# Subset FOXK2 peaks and convert rownames to GRanges
FOXK2_peaks_matrix <- motif_annotations_matrix[motif_annotations_matrix[,"FOXK2"]==TRUE,c(1:746)]

# Use custom function to extract co-ordinates as GRanges object 
FOXK2_peaks <- convert_coordinates_to_GRanges(rownames(FOXK2_peaks_matrix)) 


##############################################################################################################################
########################################## READ IN ACCESSIBILITY MODULES #####################################################

Accessibility_modules <- readLines("Supplementary_table_2_AM_peak_coordinates.txt")

# RE-FORMAT FILE INTO VECTOR LIST, WHERE EACH VECTOR IS A DIFFERENT AM
# Split each line at the semicolon
split_lines <- strsplit(Accessibility_modules, ";")

# Build the named list
AM_list <- lapply(split_lines, function(x) {
  # x[2] is the values string, split it by commas
  vals <- strsplit(x[2], ",")[[1]]
  # Clean up any extra whitespace
  vals <- trimws(vals)
  return(vals)
})

# Set the names of the list
names(AM_list) <- sapply(split_lines, function(x) trimws(x[1]))

# Get GRanges for each AM
AM1_Gr <- convert_coordinates_to_GRanges(AM_list$FullData_PM1)
AM2_Gr <- convert_coordinates_to_GRanges(AM_list$FullData_PM2)
AM3_Gr <- convert_coordinates_to_GRanges(AM_list$FullData_PM3)
AM4_Gr <- convert_coordinates_to_GRanges(AM_list$FullData_PM4)


##############################################################################################################################
################################################# CHECK OVERLAP ##############################################################
AM1_FOXK2_peaks <- AM1_Gr[AM1_Gr %in% FOXK2_peaks]
AM1_overlaps <- findOverlaps(AM1_Gr, FOXK2_peaks)
AM1_Gr[unique(queryHits(AM1_overlaps))]

AM2_FOXK2_peaks <- AM2_Gr[AM2_Gr %in% FOXK2_peaks]
AM3_FOXK2_peaks <- AM3_Gr[AM3_Gr %in% FOXK2_peaks]
AM4_FOXK2_peaks <- AM4_Gr[AM4_Gr %in% FOXK2_peaks]

summary(AM1_FOXK2_peaks) # 4/223 ranges
summary(AM2_FOXK2_peaks) # 0/63 ranges
summary(AM3_FOXK2_peaks) # 2/80 ranges
summary(AM4_FOXK2_peaks) # 8/192 ranges


# Scanning for FOXK2 motif in AMs to double check that it is largely absent

# Define background nucleotide frequencies (based on full ArchR peakset)
ACTG_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Get the FOXK2 motif from JASPAR2020
FOXK2_tfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM", name = "FOXK2"))[[1]]

# Scan for motif occurrences in the peakset
FOXK2_AM1_motif_matches <- matchMotifs(FOXK2_tfm, AM1_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                               out = "matches", p.cutoff = 1e-05, bg = ACTG_freqs)
FOXK2_AM2_motif_matches <- matchMotifs(FOXK2_tfm, AM2_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-05, bg = ACTG_freqs)
FOXK2_AM3_motif_matches <- matchMotifs(FOXK2_tfm, AM3_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-05, bg = ACTG_freqs)
FOXK2_AM4_motif_matches <- matchMotifs(FOXK2_tfm, AM4_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-05, bg = ACTG_freqs)

# Count the total hits
sum(assays(FOXK2_AM1_motif_matches)$motifMatches) # 4
sum(assays(FOXK2_AM2_motif_matches)$motifMatches) # 0
sum(assays(FOXK2_AM3_motif_matches)$motifMatches) # 2
sum(assays(FOXK2_AM4_motif_matches)$motifMatches) # 8

# Repeating the above with lower motif match threshold (1e-04)
FOXK2_AM1_motif_matches_e4 <- matchMotifs(FOXK2_tfm, AM1_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-04, bg = ACTG_freqs)
FOXK2_AM2_motif_matches_e4 <- matchMotifs(FOXK2_tfm, AM2_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-04, bg = ACTG_freqs)
FOXK2_AM3_motif_matches_e4 <- matchMotifs(FOXK2_tfm, AM3_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-04, bg = ACTG_freqs)
FOXK2_AM4_motif_matches_e4 <- matchMotifs(FOXK2_tfm, AM4_Gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                                       out = "matches", p.cutoff = 1e-04, bg = ACTG_freqs)

# Count the total hits
sum(assays(FOXK2_AM1_motif_matches_e4)$motifMatches) # 31
sum(assays(FOXK2_AM2_motif_matches_e4)$motifMatches) # 3
sum(assays(FOXK2_AM3_motif_matches_e4)$motifMatches) # 15
sum(assays(FOXK2_AM4_motif_matches_e4)$motifMatches) # 30


##############################################################################################################################
################################################# PEAKS 2 GENES ##############################################################

# Load FullData ArchR project with peak2gene links
FullData <- loadArchRProject(path = '/data/scATAC_paper/FullData_Save-ArchR')

# Get peak2GeneLinks use corCutOff = 0 to retrieve any positively correlated genes (This is what was used for GRNi. Less stringent)
p2g <- getPeak2GeneLinks(FullData, corCutOff = 0.45, FDRCutOff = 1e-04, returnLoops = FALSE)

# Convert peak2gene links dataframe to a GRanges object for further cross comparisons
p2g_granges <- metadata(p2g)$peakSet[p2g$idxATAC]
mcols(p2g_granges)$linked_gene <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]

# Getting all annotated enhancers for FOXK2 targets
FOXK2_p2G <- p2g_granges[p2g_granges %in% FOXK2_peaks]
sort(unique(FOXK2_p2G$linked_gene)) # List all predicted downstream targets of FOXK2

# Extract FOXK2 target enhancers for key placodal factors
FOXK2_targets <- c("TFAP2C", "SIX1", "DLX5", "DLX6", "GATA2", "FOXK2")  
FOXK2_target_enhancers <- FOXK2_p2G[FOXK2_p2G$linked_gene %in% FOXK2_targets]
FOXK2_targets


# Check downstream targets of different placodal AMs
AM1_targets <- p2g_granges[p2g_granges %in% AM1_Gr]
AM2_targets <- p2g_granges[p2g_granges %in% AM2_Gr]
AM3_targets <- p2g_granges[p2g_granges %in% AM3_Gr]
AM4_targets <- p2g_granges[p2g_granges %in% AM4_Gr]

write_lines(sort(unique(AM1_targets$linked_gene)), '/data/Sox8_binding_partner_analysis/scATACseq_objects/AM_target_genes/AM1_P2G_genes.txt')
write_lines(sort(unique(AM2_targets$linked_gene)), '/data/Sox8_binding_partner_analysis/scATACseq_objects/AM_target_genes/AM2_P2G_genes.txt')
write_lines(sort(unique(AM3_targets$linked_gene)), '/data/Sox8_binding_partner_analysis/scATACseq_objects/AM_target_genes/AM3_P2G_genes.txt')
write_lines(sort(unique(AM4_targets$linked_gene)), '/data/Sox8_binding_partner_analysis/scATACseq_objects/AM_target_genes/AM4_P2G_genes.txt')

# Plot a peak2gene heatmap
FullData@cellColData$condition <- paste(FullData$stage, FullData$scHelper_cell_type, sep = "_")

# Editing updated seATAC and seRNA paths to p2g metadata

S4Vectors::metadata(S4Vectors::metadata(FullData@peakSet)$Peak2GeneLinks)$seATAC <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/TransferLabel_Save-ArchR/Peak2GeneLinks/seATAC-Group-KNN.rds"
S4Vectors::metadata(S4Vectors::metadata(FullData@peakSet)$Peak2GeneLinks)$seRNA <- "/data/Sox8_binding_partner_analysis/scATACseq_objects/TransferLabel_Save-ArchR/Peak2GeneLinks/seRNA-Group-KNN.rds"

p <- plotPeak2GeneHeatmap(ArchRProj = FullData, groupBy = "stage")

readRDS("/data/scATAC_paper/FullData_Save-ArchR/Peak2GeneLinks/seATAC-Group-KNN.rds")
