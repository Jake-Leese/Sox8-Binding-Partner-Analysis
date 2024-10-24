# Script to calculate correlation coefficients and write to files

setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(TFBSTools)
library(JASPAR2024)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(universalmotif)

source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_HM_Clusters.R")

# Load and check the Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

# Subset Seurat object using custom subset_seurat_features() function that takes arguments
# Seurat_obj, DB (JASPAR2024 or AnimalTFDB) and taxa ("vertebrates" or "Homo sapiens")

HHall_ectoderm_TFs <- subset_seurat_features(HHall_ectoderm, DB = "JASPAR2024")
HHall_ectoderm_TFs_AnimalTFDB <- subset_seurat_features(HHall_ectoderm, DB = "AnimalTFDB")

#### EXTRACT COUNT DATA FOR STAGES HH8-HH9 AND CELL TYPES NC AND PLACODAL ####
