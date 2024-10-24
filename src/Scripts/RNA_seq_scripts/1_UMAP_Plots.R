# 240807 script written to produce SOX8 feature plot and UMAPs labeled by stage and ectoderm cell type for comparison

setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)

# Load and check the Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")


# view available dimensionality reductions to plot from the dataset
HHall_ectoderm@reductions

# Generate a feature plot for gene of interest
FeaturePlot(object = HHall_ectoderm, features = "SOX8", reduction = "umap.HHall", 
            cols = c("lightgrey", "red"), pt.size = 0.1)

# Plot UMAP labeled by 'stage' 
DimPlot(object = HHall_ectoderm, reduction = "umap.HHall", group.by = "stage", 
        pt.size = 0.1)

# Plot UMAP labeled by 'Ectoderm_type' 
DimPlot(object = HHall_ectoderm, reduction = "umap.HHall", group.by = "ectoderm_type", 
        pt.size = 0.1)

