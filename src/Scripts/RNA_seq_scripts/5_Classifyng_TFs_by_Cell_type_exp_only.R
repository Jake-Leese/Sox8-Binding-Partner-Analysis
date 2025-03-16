# Script to run differential gene expression analysis between TFs and identify those that are differentially expressed between PPR and NC at each stage
# In contrast to other scripts, this does not take SOX8 expression into account
# Outputs a MEME formatted motif matrix list for PPR and NC TFs at each stage
# Run using ArchR_Seurat_R_4.4.1.sif container

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
library(EnhancedVolcano)
library(DESeq2)
library(SingleCellExperiment)
library(dplyr)

source("/data/Sox8_binding_partner_analysis/Functions/Extract_HM_Clusters.R")

# Load and check the Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

head(HHall_ectoderm)

#### SUBSET THE SEURAT OBJECT INTO HH7 AND HH8 STAGES AND WITH JASPAR2024 TFs only ####

# Download JASPAR2024 motif list for vertebrates
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

# Subsetting Seurat object to include JASPAR TFs only
HHall_ectoderm_TFs <- DietSeurat(HHall_ectoderm, features = name_vector_vert)

# Extract and subset cells based on stage (HH7-HH9) 
HH7_ectoderm_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$stage %in% "HH7"]
HH8_ectoderm_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$stage %in% "HH8"]
HH9_ectoderm_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$stage %in% "HH9"]


#### PERFORM DIFFERENTIAL EXPRESSION ANALYSIS BETWEEN NC AND PLACODAL CELL TYPES ####

# Method 1: A simple Wilcox DE test that outputs a dataframe with several columns: 
# avg_log2FC: Average log2 fold change (positive for higher expression in ident.1, negative for higher expression in ident.2).
# pct.1: Proportion of cells in Placode (ident.1) expressing the gene.
# pct.2: Proportion of cells in NC (itent.2) expressing the gene.
HH7_markers_Placodal_vs_NC <- FindMarkers(
  HH7_ectoderm_TFs,
  group.by = "ectoderm_type",
  ident.1 = "placode",
  ident.2 = "NC",
  assay = "SCT",
  min.pct = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)

# Make a volcano plot of Placodal vs NC gene expression for HH7
EnhancedVolcano(HH7_markers_Placodal_vs_NC,
                lab = rownames(HH7_markers_Placodal_vs_NC),        # Gene labels (rownames)
                x = 'avg_log2FC',                      # X-axis: log2 fold change
                y = 'p_val_adj',                       # Y-axis: adjusted p-value
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "adjusted p-value"),
                pCutoff = 0.05,                        # Significance threshold
                FCcutoff = 0.25,                       # Fold change threshold
                title = 'HH7 Volcano Plot: Placodal vs NC',
                labSize = 3.0,                         # Label size
                pointSize = 3.0,                       # Point size
                col = c('black', 'blue', 'green', 'red'))  # Colors for points


# Extract genes expressed significantly more in Placode or NC cells, and those with relatively high expression in both. + log2FC is placodal, - log2FC is NC
HH7_Placodal_TFs <- rownames(HH7_markers_Placodal_vs_NC[HH7_markers_Placodal_vs_NC$avg_log2FC > 0 & HH7_markers_Placodal_vs_NC$p_val_adj < 0.05, ])
HH7_NC_TFs <- rownames(HH7_markers_Placodal_vs_NC[HH7_markers_Placodal_vs_NC$avg_log2FC < 0 & HH7_markers_Placodal_vs_NC$p_val_adj < 0.05, ])
HH7_NC_and_Placodal_TFs <- rownames(HH7_markers_Placodal_vs_NC[(HH7_markers_Placodal_vs_NC$pct.1 > 0.2 & HH7_markers_Placodal_vs_NC$pct.2 > 0.2 & HH7_markers_Placodal_vs_NC$p_val_adj > 0.05) |
                                                               (HH7_markers_Placodal_vs_NC$pct.1 > 0.4 & HH7_markers_Placodal_vs_NC$pct.2 > 0.4), ])
                                                               

# Subset the DE test results based on TF lists. Write to file
HH7_Placodal_TFs <- HH7_markers_Placodal_vs_NC[rownames(HH7_markers_Placodal_vs_NC) %in% HH7_Placodal_TFs, ]
HH7_NC_TFs <- HH7_markers_Placodal_vs_NC[rownames(HH7_markers_Placodal_vs_NC) %in% HH7_NC_TFs, ]
HH7_NC_and_Placodal_TFs <- HH7_markers_Placodal_vs_NC[rownames(HH7_markers_Placodal_vs_NC) %in% HH7_NC_and_Placodal_TFs, ]

write.csv(HH7_Placodal_TFs, "DE_analyses/HH7_Placodal_TFs_wilcox.csv")
write.csv(HH7_NC_TFs, "DE_analyses/HH7_NC_TFs_wilcox.csv")
write.csv(HH7_NC_and_Placodal_TFs, "DE_analyses/HH7_NC_and_Placodal_TFs_wilcox.csv")


# Repeat above steps for HH8 stage
HH8_markers_Placodal_vs_NC <- FindMarkers(
  HH8_ectoderm_TFs,
  group.by = "ectoderm_type",
  ident.1 = "placode",
  ident.2 = "NC",
  assay = "SCT",
  min.pct = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)

# Make a volcano plot of Placodal vs NC gene expression for HH8
EnhancedVolcano(HH8_markers_Placodal_vs_NC,
                lab = rownames(HH8_markers_Placodal_vs_NC),        # Gene labels (rownames)
                x = 'avg_log2FC',                      # X-axis: log2 fold change
                y = 'p_val_adj',                       # Y-axis: adjusted p-value
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "adjusted p-value"),
                pCutoff = 0.05,                        # Significance threshold
                FCcutoff = 0.25,                       # Fold change threshold
                title = 'HH8 Volcano Plot: Placodal vs NC',
                labSize = 3.0,                         # Label size
                pointSize = 3.0,                       # Point size
                col = c('black', 'blue', 'green', 'red'))  # Colors for points


# Extract genes expressed significantly more in Placode or NC cells, and those with relatively high expression in both. + log2FC is placodal, - log2FC is NC
HH8_Placodal_TFs <- rownames(HH8_markers_Placodal_vs_NC[HH8_markers_Placodal_vs_NC$avg_log2FC > 0 & HH8_markers_Placodal_vs_NC$p_val_adj < 0.05, ])
HH8_NC_TFs <- rownames(HH8_markers_Placodal_vs_NC[HH8_markers_Placodal_vs_NC$avg_log2FC < 0 & HH8_markers_Placodal_vs_NC$p_val_adj < 0.05, ])
HH8_NC_and_Placodal_TFs <- rownames(HH8_markers_Placodal_vs_NC[(HH8_markers_Placodal_vs_NC$pct.1 > 0.2 & HH8_markers_Placodal_vs_NC$pct.2 > 0.2 & HH8_markers_Placodal_vs_NC$p_val_adj > 0.05) |
                                                                 (HH8_markers_Placodal_vs_NC$pct.1 > 0.4 & HH8_markers_Placodal_vs_NC$pct.2 > 0.4), ])


# Subset the DE test results based on TF lists. Write to file
HH8_Placodal_TFs <- HH8_markers_Placodal_vs_NC[rownames(HH8_markers_Placodal_vs_NC) %in% HH8_Placodal_TFs, ]
HH8_NC_TFs <- HH8_markers_Placodal_vs_NC[rownames(HH8_markers_Placodal_vs_NC) %in% HH8_NC_TFs, ]
HH8_NC_and_Placodal_TFs <- HH8_markers_Placodal_vs_NC[rownames(HH8_markers_Placodal_vs_NC) %in% HH8_NC_and_Placodal_TFs, ]

write.csv(HH8_Placodal_TFs, "DE_analyses/HH8_Placodal_TFs_wilcox.csv")
write.csv(HH8_NC_TFs, "DE_analyses/HH8_NC_TFs_wilcox.csv")
write.csv(HH8_NC_and_Placodal_TFs, "DE_analyses/HH8_NC_and_Placodal_TFs_wilcox.csv")


# Repeat above steps for HH9 stage
HH9_markers_Placodal_vs_NC <- FindMarkers(
  HH9_ectoderm_TFs,
  group.by = "ectoderm_type",
  ident.1 = "placode",
  ident.2 = "NC",
  assay = "SCT",
  min.pct = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)

# Make a volcano plot of Placodal vs NC gene expression for HH9
EnhancedVolcano(HH9_markers_Placodal_vs_NC,
                lab = rownames(HH9_markers_Placodal_vs_NC),        # Gene labels (rownames)
                x = 'avg_log2FC',                      # X-axis: log2 fold change
                y = 'p_val_adj',                       # Y-axis: adjusted p-value
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "adjusted p-value"),
                pCutoff = 0.05,                        # Significance threshold
                FCcutoff = 0.25,                       # Fold change threshold
                title = 'HH9 Volcano Plot: Placodal vs NC',
                labSize = 3.0,                         # Label size
                pointSize = 3.0,                       # Point size
                col = c('black', 'blue', 'green', 'red'))  # Colors for points


# Extract genes expressed significantly more in Placode or NC cells, and those with relatively high expression in both. + log2FC is placodal, - log2FC is NC
HH9_Placodal_TFs <- rownames(HH9_markers_Placodal_vs_NC[HH9_markers_Placodal_vs_NC$avg_log2FC > 0 & HH9_markers_Placodal_vs_NC$p_val_adj < 0.05, ])
HH9_NC_TFs <- rownames(HH9_markers_Placodal_vs_NC[HH9_markers_Placodal_vs_NC$avg_log2FC < 0 & HH9_markers_Placodal_vs_NC$p_val_adj < 0.05, ])
HH9_NC_and_Placodal_TFs <- rownames(HH9_markers_Placodal_vs_NC[(HH9_markers_Placodal_vs_NC$pct.1 > 0.2 & HH9_markers_Placodal_vs_NC$pct.2 > 0.2 & HH9_markers_Placodal_vs_NC$p_val_adj > 0.05) |
                                                                 (HH9_markers_Placodal_vs_NC$pct.1 > 0.4 & HH9_markers_Placodal_vs_NC$pct.2 > 0.4), ])


# Subset the DE test results based on TF lists. Write to file
HH9_Placodal_TFs <- HH9_markers_Placodal_vs_NC[rownames(HH9_markers_Placodal_vs_NC) %in% HH9_Placodal_TFs, ]
HH9_NC_TFs <- HH9_markers_Placodal_vs_NC[rownames(HH9_markers_Placodal_vs_NC) %in% HH9_NC_TFs, ]
HH9_NC_and_Placodal_TFs <- HH9_markers_Placodal_vs_NC[rownames(HH9_markers_Placodal_vs_NC) %in% HH9_NC_and_Placodal_TFs, ]

write.csv(HH9_Placodal_TFs, "DE_analyses/HH9_Placodal_TFs_wilcox.csv")
write.csv(HH9_NC_TFs, "DE_analyses/HH9_NC_TFs_wilcox.csv")
write.csv(HH9_NC_and_Placodal_TFs, "DE_analyses/HH9_NC_and_Placodal_TFs_wilcox.csv")



#### PREPARE MEME FORMAT MOTIF LISTS BASED OFF CO-EXPRESSION ANALYSIS ABOVE ####

# Read in SOX8 co-expression heatmap clusters and extract TFs based on their co-expression patterns
HH7_Placodal_TFs <- read.csv("DE_analyses/HH7_Placodal_TFs_wilcox.csv")
HH7_Placodal_TFs <- HH7_Placodal_TFs$X
HH7_NC_TFs <- read.csv("DE_analyses/HH7_NC_TFs_wilcox.csv")
HH7_NC_TFs <- HH7_NC_TFs$X
HH7_NC_and_Placodal_TFs <- read.csv("DE_analyses/HH7_NC_and_Placodal_TFs_wilcox.csv")
HH7_NC_and_Placodal_TFs <- HH7_NC_and_Placodal_TFs$X
HH8_Placodal_TFs <- read.csv("DE_analyses/HH8_Placodal_TFs_wilcox.csv")
HH8_Placodal_TFs <- HH8_Placodal_TFs$X
HH8_NC_TFs <- read.csv("DE_analyses/HH8_NC_TFs_wilcox.csv")
HH8_NC_TFs <- HH8_NC_TFs$X
HH8_NC_and_Placodal_TFs <- read.csv("DE_analyses/HH8_NC_and_Placodal_TFs_wilcox.csv")
HH8_NC_and_Placodal_TFs <- HH8_NC_and_Placodal_TFs$X
HH9_Placodal_TFs <- read.csv("DE_analyses/HH9_Placodal_TFs_wilcox.csv")
HH9_Placodal_TFs <- HH9_Placodal_TFs$X
HH9_NC_TFs <- read.csv("DE_analyses/HH9_NC_TFs_wilcox.csv")
HH9_NC_TFs <- HH9_NC_TFs$X
HH9_NC_and_Placodal_TFs <- read.csv("DE_analyses/HH9_NC_and_Placodal_TFs_wilcox.csv")
HH9_NC_and_Placodal_TFs <- HH9_NC_and_Placodal_TFs$X

# subset motif list based on motif lists of interest
HH7_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_Placodal_TFs)]
HH7_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_NC_TFs)]
HH7_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH7_NC_and_Placodal_TFs)]
HH8_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_Placodal_TFs)]
HH8_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_NC_TFs)]
HH8_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH8_NC_and_Placodal_TFs)]
HH9_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH9_Placodal_TFs)]
HH9_NC_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH9_NC_TFs)]
HH9_NC_and_Placodal_MotifList_vert <- motifList_vert[intersect(names(motifList_vert), HH9_NC_and_Placodal_TFs)]

# Extract SOX8 motif matrix to use as primary motif in SpaMo analysis
Sox8_motif <- motifList_vert$SOX8

# Full scATAC peakset background ACTG_freqs = 0.2487819, 0.2512324, 0.2512223, 0.2487634
ACGT_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Write JASPAR Motif matrices to minimal meme format

write_meme(HH7_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH7_Placodal_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH7_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH7_NC_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH7_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH7_NC_and_Placodal_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH8_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_Placodal_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH8_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_NC_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH8_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_NC_and_Placodal_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH9_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH9_Placodal_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH9_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH9_NC_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH9_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH9_NC_and_Placodal_motifs_wilcox_PWM.txt", version = 5, ACGT_freqs)




# Method 2 for differential expression analysis. Make a pseudobulk dataset and use DESeq2 #
# Note that, for the sake of simply establishing whether genes are expressed in one cell type or the other, or both, this method is not
# suitable due to the large difference in cell numnbers between NC and placodal cells. DESeq2 requires aggregating the counts, and by this
# method will greatly bias placode cell type expression, due to the ~8 fold higher number of cells. Leaving code here in case I come back to it


# First extract SCT count data and cell metadata from Seurat object and filter to NC and Placode cells only
HH7_Metadata <- HH7_ectoderm_TFs@meta.data  # 4036 cells
HH7_Metadata_filtered <- HH7_Metadata %>% filter(ectoderm_type %in% c("NC", "placode")) # 2589 cells

HH7_SCT_counts <- GetAssayData(HH7_ectoderm_TFs, layer = "counts", assay = "SCT")
HH7_SCT_counts_filtered <- HH7_SCT_counts[, colnames(HH7_SCT_counts) %in% rownames(HH7_Metadata_filtered)]

# Create a pseudo-bulk matrix by calculating sum counts for each cell type
HH7_pseudo_bulk <- as.data.frame(t(HH7_SCT_counts_filtered)) %>%
  mutate(cell_type = HH7_Metadata_filtered$ectoderm_type) %>%
  group_by(cell_type) %>%
  summarise(across(everything(), sum))

# Scale the counts by the number of cells in each group (proportional counts) (NC cells = 306, Placode cells = 2283)
nc_cells <- nrow(HH7_Metadata_filtered[HH7_Metadata_filtered$ectoderm_type == "NC", ])
placode_cells <- nrow(HH7_Metadata_filtered[HH7_Metadata_filtered$ectoderm_type == "placode", ])

pseudo_bulk_scaled <- HH7_pseudo_bulk %>%
  mutate(across(starts_with("NC"), ~ . / nc_cells),
         across(starts_with("placode"), ~ . / placode_cells)) # Didn't do anything?

# Transpose to get genes as rows and pseudo-bulk samples as columns
HH7_pseudo_bulk <- t(as.matrix(HH7_pseudo_bulk[, -1]))
colnames(HH7_pseudo_bulk) <- c("NC", "Placode")

# Create a data frame for the DESeq2 design (experimental conditions)
condition <- factor(c("NC", "Placode"))

# Construct DESeq2 dataset
HH7_DESeq2_obj <- DESeqDataSetFromMatrix(countData = HH7_pseudo_bulk,
                              colData = data.frame(condition = condition),
                              design = ~ condition)
