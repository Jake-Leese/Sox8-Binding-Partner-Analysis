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

# Pull PWM motif matrix from JASPAR2020
motifList_vert <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))
name_vector_vert <- c()
for (i in 1:length(motifList_vert)){
  name <- name(motifList_vert[[i]])
  name_vector_vert <- c(name_vector_vert, name)
}
names(motifList_vert) <- name_vector_vert

# annotate peaks in ArchR object with these motifs
ss4ArchRProj <- addMotifAnnotations(ss4ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)
ss8ArchRProj <- addMotifAnnotations(ss8ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)

# Pull motif positions from peaksets
ss4_motifPositions <- getPositions(ss4ArchRProj)
ss8_motifPositions <- getPositions(ss8ArchRProj)

motifs <- c("SOX8", "LMX1A", "PRDM1", "DLX6")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(ss4_motifPositions), value = TRUE)))

# Adding broad cell type group coverages
ss4ArchRProj <- addGroupCoverages(ss4ArchRProj, groupBy = "transferred_scHelper_cell_type_broad")
ss8ArchRProj <- addGroupCoverages(ss8ArchRProj, groupBy = "transferred_scHelper_cell_type_broad")

# Footprinting markerMotifs
ss4Footprints <- getFootprints(
  ss4ArchRProj,
  positions = ss4_motifPositions[markerMotifs],
  groupBy = "transferred_scHelper_cell_type_broad"
)

ss8Footprints <- getFootprints(
  ss8ArchRProj,
  positions = ss8_motifPositions[markerMotifs],
  groupBy = "transferred_scHelper_cell_type_broad"
)

plotFootprints(
  seFoot = ss4Footprints,
  ArchRProj = ss4ArchRProj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = ss8Footprints,
  ArchRProj = ss8ArchRProj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
