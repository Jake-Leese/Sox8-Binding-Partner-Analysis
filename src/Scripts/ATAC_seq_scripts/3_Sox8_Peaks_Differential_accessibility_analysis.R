# Script to perform differential accessibility analysis, make volcano plots and extract differentially accessible peaks

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

# Load custom R functions
source("/data/Sox8_binding_partner_analysis/Functions/ArchRAddUniqueIdsToSe.R")
source("/data/Sox8_binding_partner_analysis/Functions/ArchR_ExtractIds.R")
source("/data/Sox8_binding_partner_analysis/Functions/AddArchRMetaData.R")

# Set number of threads
addArchRThreads(threads = 64)

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Ggallus.UCSC.galGal6)

geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Ggallus.UCSC.galGal6.refGene, 
                                       OrgDb = org.Gg.eg.db)

# Loading ArchR projects
ss4_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR/', force = FALSE)
ss8_Sox8_ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Sox8_Save_ArchR/', force = FALSE)


#### PERFORM DIFFERENTIAL ACCESSIBILITY ANALYSES AND VISUALIZE RESULTS IN HEATMAPS AND VOLCANO PLOTS ####

# Create a summarized experiment to identify marker peaks for different cell groups, in this case broad ectodermal classifications.
ss8_markerPeaks <- getMarkerFeatures(
  ArchRProj = ss8_Sox8_ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

# Plotting heatmap of cluster peaks
heatmapPeaks <- markerHeatmap(
  seMarker = ss8_markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# A pairwise differential accessibility analysis between Placodal and NC cell types
ss8_markerTest <- getMarkerFeatures(
  ArchRProj = ss8_Sox8_ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "NC")

# This function transfers metadata from ArchR project matrix (PeakMatrix in this case),
# and it assigns the corresponding peak in the getMarkerFeatures() object. 
# Specifically it pulls over the "nearestGene" and unique_id as "seqnames":"start"-"end"
ss8_markerTest <- ArchRAddUniqueIdsToSe(ss8_markerTest, ss8_Sox8_ArchRProj,"PeakMatrix")

Peak_IDs <- ArchR_ExtractIds(ss8_markerTest, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1")

# Plotting a volcano plot between Sox8 target accessibility in placodal vs NC cell types
SS8_pv <- plotMarkers(seMarker = ss8_markerTest, 
                  name = "Placodal", 
                  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", 
                  plotAs = "Volcano")

plot(SS8_pv)

# Save .svg version of this plot
svg("ss8_NC_v_Placode_volcano.svg", width = 8, height = 6)
plot(pv)
dev.off()

# Save vectorized version of this plot. Up-regulated = Placodal, Down-regulated = NC.
plotPDF(pv, 
        name = "Placodal-vs-NC-SOX8-positive-Markers-Volcano", 
        width = 5, 
        height = 5, 
        ArchRProj = ss8_Sox8_ArchRProj, 
        addDOC = FALSE)

# PLOTTING SS4 DIFFERENTIAL ACCESSIBILITY VOLCANO PLOT
ss4_markerTest <- getMarkerFeatures(
  ArchRProj = ss4_Sox8_ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "NC")

SS4_pv <- plotMarkers(seMarker = ss4_markerTest, 
                      name = "Placodal", 
                      cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", 
                      plotAs = "Volcano")

plot(SS4_pv)

svg("ss4_NC_v_Placode_volcano.svg", width = 8, height = 6)
plot(SS4_pv)
dev.off()

#### EXTRACT PEAKS BASED OF ACCESSIBILITY ####

# Extract peaks with higher accessibility in Placodal cells
ss8_Sox8_Placodal_peaks <- getMarkers(ss8_markerTest,
           cutOff = "FDR <= 0.1 & (Log2FC) >= 1",
           returnGR = TRUE)
ss8_Sox8_Placodal_peaks <- ss8_Sox8_Placodal_peaks$Placodal # Extracts peaks as GRanges object from list

# Extract peaks with higher accessibility in Neural Crest cells
ss8_Sox8_NC_peaks <- getMarkers(ss8_markerTest,
                      cutOff = "FDR <= 0.1 & (Log2FC) <= -1",
                      returnGR = TRUE)
ss8_Sox8_NC_peaks <- ss8_Sox8_NC_peaks$Placodal

# Extract peaks with high accessibility over both cell types
ss8_Sox8_NC_and_Placodal_peaks <- getMarkers(ss8_markerTest,
  cutOff = "FDR > 0.1 & abs(Log2FC) < 1",  
  returnGR = TRUE)

ss8_Sox8_NC_and_Placodal_peaks <- ss8_Sox8_NC_and_Placodal_peaks$Placodal

# add nearest genes to GRanges objects as additional metadata column
ss8_Sox8_Placodal_peaks <- AddArchRMetaData(ss8_Sox8_Placodal_peaks, ss8_Sox8_ArchRProj)
ss8_Sox8_NC_peaks <- AddArchRMetaData(ss8_Sox8_NC_peaks, ss8_Sox8_ArchRProj)
ss8_Sox8_NC_and_Placodal_peaks <- AddArchRMetaData(ss8_Sox8_NC_and_Placodal_peaks, ss8_Sox8_ArchRProj)

# Save GRanges objects as rds files
saveRDS(ss8_Sox8_Placodal_peaks, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_Placodal_peaks.rds")
saveRDS(ss8_Sox8_NC_peaks, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_peaks.rds")
saveRDS(ss8_Sox8_NC_and_Placodal_peaks, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss8_Sox8_NC_and_Placodal_peaks.rds")


#### REPEATING ABOVE FOR SS4 STAGE ####

# Plotting a heatmap of peak accessibility by broad cell type labels
ss4_markerPeaks <- getMarkerFeatures(
  ArchRProj = ss4_Sox8_ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

# Plotting heatmap of cluster peaks
ss4_heatmapPeaks <- plotMarkerHeatmap(
  seMarker = ss4_markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE)

draw(ss4_heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# A pairwise differential accessibility analysis between Placodal and NC cell types
ss4_markerTest <- getMarkerFeatures(
  ArchRProj = ss4_Sox8_ArchRProj,
  useMatrix = "PeakMatrix",
  groupBy = "transferred_scHelper_cell_type_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Placodal",
  bgdGroups = "NC")

# This function transfers metadata from ArchR project matrix (PeakMatrix in this case),
# and it assigns the corresponding peak in the getMarkerFeatures() object. 
# Specifically it pulls over the "nearestGene" and unique_id as "seqnames":"start"-"end"
ss4_markerTest <- ArchRAddUniqueIdsToSe(ss4_markerTest, ss4_Sox8_ArchRProj,"PeakMatrix")

Peak_IDs <- ArchR_ExtractIds(ss4_markerTest, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1")

# Plotting a volcano plot between Sox8 target accessibility in placodal vs NC cell types
pv <- plotMarkers(seMarker = ss4_markerTest, 
                  name = "Placodal", 
                  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", 
                  plotAs = "Volcano")

plot(pv)


# Extract peaks based on accessibility across placodal and NC cells
ss4_Sox8_Placodal_peaks <- getMarkers(ss4_markerTest,
                                      cutOff = "FDR <= 0.1 & (Log2FC) >= 1",
                                      returnGR = TRUE)
ss4_Sox8_Placodal_peaks <- ss4_Sox8_Placodal_peaks$Placodal # Extracts peaks as GRanges object from list


ss4_Sox8_NC_peaks <- getMarkers(ss4_markerTest,
                                cutOff = "FDR <= 0.1 & (Log2FC) <= -1",
                                returnGR = TRUE)
ss4_Sox8_NC_peaks <- ss4_Sox8_NC_peaks$Placodal


ss4_Sox8_NC_and_Placodal_peaks <- getMarkers(ss4_markerTest,
                                             cutOff = "FDR > 0.1 & abs(Log2FC) < 1",  
                                             returnGR = TRUE)

ss4_Sox8_NC_and_Placodal_peaks <- ss4_Sox8_NC_and_Placodal_peaks$Placodal

# add ArchR peakset metadata to GRanges
ss4_Sox8_Placodal_peaks <- AddArchRMetaData(ss4_Sox8_Placodal_peaks, ss4_Sox8_ArchRProj)
ss4_Sox8_NC_peaks <- AddArchRMetaData(ss4_Sox8_NC_peaks, ss4_Sox8_ArchRProj)
ss4_Sox8_NC_and_Placodal_peaks <- AddArchRMetaData(ss4_Sox8_NC_and_Placodal_peaks, ss4_Sox8_ArchRProj)

# Check that there is no overlap between any of the GRanges. They should all be mutually exclusive
findOverlaps(ss4_Sox8_Placodal_peaks, ss4_Sox8_NC_and_Placodal_peaks) # 0 HITS
findOverlaps(ss4_Sox8_Placodal_peaks, ss4_Sox8_NC_peaks)              # 0 HITS
findOverlaps(ss4_Sox8_NC_and_Placodal_peaks, ss4_Sox8_NC_peaks)       # 0 HITS

# Save GRanges objects as rds files
saveRDS(ss4_Sox8_Placodal_peaks, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_Placodal_peaks.rds")
saveRDS(ss4_Sox8_NC_peaks, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_peaks.rds")
saveRDS(ss4_Sox8_NC_and_Placodal_peaks, "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/ss4_Sox8_NC_and_Placodal_peaks.rds")





#### MOTIF ENRICHMENT ANALYSIS ####

# Checking for motifs enriched in up-regulated peaks (placodal accessible)
Placodal_motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = ss8_Sox8_ArchRProj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# Prepare summarized experiment for plotting with ggplot
df <- data.frame(TF = rownames(Placodal_motifsUp), mlog10Padj = assay(Placodal_motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

# plotting
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
