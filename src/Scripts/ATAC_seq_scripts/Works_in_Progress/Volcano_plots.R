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


# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)


# Perform differential accessibility analysis
ss8_PPR_v_NC_markers <- getMarkerFeatures(
  ArchRProj = ss8ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "NC"
)

ss8_PPR_v_Neural_markers <- getMarkerFeatures(
  ArchRProj = ss8ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "Neural"
)

ss8_PPR_v_Contam_markers <- getMarkerFeatures(
  ArchRProj = ss8ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "Contam"
)

ss8_PPR_marker_peaks <- getMarkerFeatures(
  ArchRProj = ss8ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Plot differential accessibility volcano plots for each pairwise comparison of HH9 PPR
ss8_PPR_marker_peaks_pv <- plotMarkers(seMarker = ss8_PPR_marker_peaks, 
                               name = "Placodal", 
                               cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", 
                               plotAs = "Volcano")

plot(ss8_PPR_marker_peaks_pv)

Plot_path <- "/data/Sox8_binding_partner_analysis/Plots/March_25_Volcano_plots/"

svg(paste0(Plot_path, "ss8_Placode_v_All_volcano.svg"), width = 8, height = 6)
plot(ss8_PPR_marker_peaks_pv)
dev.off()
