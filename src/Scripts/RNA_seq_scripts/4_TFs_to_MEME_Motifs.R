# Script to prepare MEME compatible motif lists for co-expressed TFs
# Using ArchR_Seurat_R_4.4.1.sif container


.libPaths("/R/libs/ArchR_Seurat_R_441")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(BSgenome.Ggallus.UCSC.galGal6)
library(JASPAR2024)
library(universalmotif)
library(TFBSTools)
library(dplyr)
library(tidyr)
library(purrr)


##### READ IN CO-EXPRESSED TFs AND PREPARE MEME COMPATIBLE MOTIF LISTS #####

# All TFs that have co-expression values with Sox8 that exceed the mean+stdev for at leasst one stage and cell type combination
All_Sox8_CoExpressed_TFs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/Heatmap_clusters/Sox8_coexpression_clusters_HH8_HH9_JASPAR_RNA_assay_mean+stdv_12C_SOX8_cells.csv")

# Subsetting TFs based on whether they have correlation scores over the mean + stdev values calculated for each cell type at each stage
All_HH8_Placodal_TFs <- All_Sox8_CoExpressed_TFs[All_Sox8_CoExpressed_TFs$HH8_Placode >= 0.109216,]$X  
All_HH8_NC_TFs <- All_Sox8_CoExpressed_TFs[All_Sox8_CoExpressed_TFs$HH8_NC >= 0.117146,]$X  
All_HH9_Placodal_TFs <- All_Sox8_CoExpressed_TFs[All_Sox8_CoExpressed_TFs$HH9_Placode >= 0.148648,]$X  
All_HH9_NC_TFs <- All_Sox8_CoExpressed_TFs[All_Sox8_CoExpressed_TFs$HH9_NC >= 0.208392,]$X

# TFs that are differentially expressed between NC and PPR at each stage
HH8_Placodal_DEGs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH8_Placodal_TFs_wilcox.csv")$X
HH8_NC_DEGs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH8_NC_TFs_wilcox.csv")$X
HH9_Placodal_DEGs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH9_Placodal_TFs_wilcox.csv")$X
HH9_NC_DEGs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH9_NC_TFs_wilcox.csv")$X

# Add SOX8 to all TF vectors so that it is included as a potential binding partner in the spamo analysis
add_SOX8 <- function(vec) {
  if (!"SOX8" %in% vec) {
    vec <- c(vec, "SOX8")
  }
  return(vec)
} # Function to check if vector has "SOX8", and if not, add it to the vector
All_HH8_Placodal_TFs <- add_SOX8(All_HH8_Placodal_TFs)
All_HH8_NC_TFs <- add_SOX8(All_HH8_NC_TFs)
All_HH9_Placodal_TFs <- add_SOX8(All_HH9_Placodal_TFs)
All_HH9_NC_TFs <- add_SOX8(All_HH9_NC_TFs)
HH8_Placodal_DEGs <- add_SOX8(HH8_Placodal_DEGs)
HH8_NC_DEGs <- add_SOX8(HH8_NC_DEGs)
HH9_Placodal_DEGs <- add_SOX8(HH9_Placodal_DEGs)
HH9_NC_DEGs <- add_SOX8(HH9_NC_DEGs)


# Prepare motif lists
Jaspar2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(Jaspar2024))
motifList_vert <- getMatrixSet(x = sq24, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PFM"))

# rename each motif from their ID to their TF name
name_vector_vert <- c()
for (i in 1:length(motifList_vert)){
  name <- name(motifList_vert[[i]])
  name_vector_vert <- c(name_vector_vert, name)
}
names(motifList_vert) <- name_vector_vert

# subset motif list based on motif lists of interest
All_HH8_Placodal_motifs <- motifList_vert[intersect(names(motifList_vert), All_HH8_Placodal_TFs)]
All_HH8_NC_motifs <- motifList_vert[intersect(names(motifList_vert), All_HH8_NC_TFs)]
All_HH9_Placodal_motifs <- motifList_vert[intersect(names(motifList_vert), All_HH9_Placodal_TFs)]
All_HH9_NC_motifs <- motifList_vert[intersect(names(motifList_vert), All_HH9_NC_TFs)]

HH8_Placodal_DEGs_motifs <- motifList_vert[intersect(names(motifList_vert), HH8_Placodal_DEGs)]
HH8_NC_DEGs_motifs <- motifList_vert[intersect(names(motifList_vert), HH8_NC_DEGs)]
HH9_Placodal_DEGs_motifs <- motifList_vert[intersect(names(motifList_vert), HH9_Placodal_DEGs)]
HH9_NC_DEGs_motifs <- motifList_vert[intersect(names(motifList_vert), HH9_NC_DEGs)]

# Extract SOX8 motif matrix to use as primary motif in SpaMo analysis
Sox8_motif <- motifList_vert$SOX8

# Full scATAC peakset background ACTG_freqs = 0.2487819, 0.2512324, 0.2512223, 0.2487634
ACGT_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)


# Write JASPAR Motif matrices to minimal meme format
write_meme(Sox8_motif, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/Sox8.txt", version = 5, ACGT_freqs, overwrite = TRUE)

write_meme(All_HH8_Placodal_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/All_HH8_Placodal_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)
write_meme(All_HH8_NC_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/All_HH8_NC_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)
write_meme(All_HH9_Placodal_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/All_HH9_Placodal_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)
write_meme(All_HH9_NC_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/All_HH9_NC_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)

write_meme(HH8_Placodal_DEGs_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/HH8_Placodal_DEGs_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)
write_meme(HH8_NC_DEGs_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/HH8_NC_DEGs_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)
write_meme(HH9_Placodal_DEGs_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/HH9_Placodal_DEGs_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)
write_meme(HH9_NC_DEGs_motifs, "/data/Sox8_binding_partner_analysis/meme_suite/March_2025/MEME_motifs/HH9_NC_DEGs_motifs.txt", version = 5, ACGT_freqs, overwrite = TRUE)




