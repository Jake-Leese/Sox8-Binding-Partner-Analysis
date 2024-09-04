# Script to setup ArchR, save/load projects, and to subset an archR project peakset, in this case based on the presence of a Sox8 binding motif

.libPaths("/R/libs/ArchR_Seurat_R_441")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(ArchR)
library(igraph)
library(Seurat)
library(TFBSTools)
library(dplyr)
library(tidyr)
library(purrr)
library(JASPAR2024)
library(BSgenome.Ggallus.UCSC.galGal6)
library(TxDb.Ggallus.UCSC.galGal6.refGene)
library(org.Gg.eg.db)

# Set number of threads
addArchRThreads(threads = 64)

## Creating genome and gene annotations for ArchR to interact with

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Ggallus.UCSC.galGal6)

geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Ggallus.UCSC.galGal6.refGene, 
                                       OrgDb = org.Gg.eg.db)

# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

# Testing ArchR project is functional by plotting a simple UMAP
p1 <- plotEmbedding(
  ArchRProj = ss4ArchRProj, 
  colorBy = "cellColData",
  name = "transferred_scHelper_cell_type_broad",
  embedding = "UMAP")

getAvailableMatrices(ss8ArchRProj)

# download motif database
Jaspar2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(Jaspar2024))
motifList_vert <- getMatrixSet(x = sq24, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))
motifList_Hs <- getMatrixSet(x = sq24, opts = list(collection = "CORE", species = "Homo sapiens", matrixtype = "PWM"))

# rename each motif from their ID to their TF name
name_vector_vert <- c()
for (i in 1:length(motifList_vert)){
  name <- name(motifList_vert[[i]])
  name_vector_vert <- c(name_vector_vert, name)
}
names(motifList_vert) <- name_vector_vert

name_vector_Hs <- c()
for (i in 1:length(motifList_Hs)){
  name <- name(motifList_Hs[[i]])
  name_vector_Hs <- c(name_vector_Hs, name)
}
names(motifList_Hs) <- name_vector_Hs

# annotate peaks in ArchR object with these motifs
ss4ArchRProj <- addMotifAnnotations(ss4ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)
ss8ArchRProj <- addMotifAnnotations(ss8ArchRProj, name = "Motif", motifPWMs = motifList_vert, cutOff = 1e-05, force = T)
print("Motifs matrix added to ArchR object!")


### EXTRACT THE MOTIF_MATCH/PEAK MATRIX FROM THE ARCHR PROJECT ###

# see where the motif x peak matrix is stored here. Matches are stored in an .rds file whose path is given in the "Matches" variable within the peak annotation object
ss4_vert_anno <- getPeakAnnotation(ss4ArchRProj, name = "Motif")
ss4_vert_anno$Matches

# need to read in this rds file to extract peak annotations with motifs - dont know if this will work in nextflow...
# motifs_peaks <- readRDS("~/output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR/Annotations/Motif-Matches-In-Peaks.rds")
ss4_vert_motif_peaks <- readRDS(ss4_vert_anno$Matches)
dim(ss4_vert_motif_peaks) 

# then need to extract the sparse matrix from the ranged experiment object
ss4_vert_motif_matrix <- assays(ss4_vert_motif_peaks)$matches
rownames(ss4_vert_motif_matrix) <- rowData(ss4_vert_motif_peaks)$name
colnames(ss4_vert_motif_matrix) <- colnames(ss4_vert_motif_peaks)

# Save the Motif X Peak matrix for future use
saveRDS(ss4_vert_motif_matrix, file = "MotifMatrices/SS4_vert_motif_matrix_1e-05.rds")


# Repeat all the above for the ss8 object
ss8_vert_anno <- getPeakAnnotation(ss8ArchRProj, name = "Motif")
ss8_vert_anno$Matches

ss8_vert_motif_peaks <- readRDS(ss8_vert_anno$Matches)
dim(ss8_vert_motif_peaks) 

ss8_vert_motif_matrix <- assays(ss8_vert_motif_peaks)$matches
rownames(ss8_vert_motif_matrix) <- rowData(ss8_vert_motif_peaks)$name
colnames(ss8_vert_motif_matrix) <- colnames(ss8_vert_motif_peaks)

saveRDS(ss8_vert_motif_matrix, file = "MotifMatrices/SS8_vert_motif_matrix_1e-05.rds")

# Saving the Arch R projects
saveArchRProject(ss4ArchRProj, "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR")
saveArchRProject(ss8ArchRProj, "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR")



