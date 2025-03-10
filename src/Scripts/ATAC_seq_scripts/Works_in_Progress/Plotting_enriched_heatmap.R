# Attempt to plot my own heatmaps Using the subsetted ss4 and ss8 peaks.
# The statistics involved to do this well is a lot more complicated than I foresaw and I am therefore parking this to one side for now

.libPaths("/R/libs/AT_ArcR_macs2")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

# Load necessary libraries
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
library(JASPAR2020)
library(motifmatchr)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)

# Define Paths to ArchR Projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)


################################################################################################################################
#################################### LOADING GRANGES FOR ALL PPR AND NC PEAKS ##################################################

# Read in full peaksets
ss4_PPR_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss4_PPR_all_Peaks.gr.rds")
ss4_NC_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss4_NC_all_Peaks.gr.rds")
ss8_PPR_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss8_PPR_all_Peaks.gr.rds")
ss8_NC_all_Peaks <- readRDS("/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/250225_updated_peaksets/ss8_NC_all_Peaks.gr.rds")

# Combine the GRanges objects and remove duplicates
NC_PPR_combined_gr <- unique(c(ss4_PPR_all_Peaks, ss4_NC_all_Peaks, ss8_PPR_all_Peaks, ss8_NC_all_Peaks))


#########################################################################################################
######################################## MOTIF SCANNING #################################################

# Define background nucleotide frequencies (based on full ArchR peakset)
ACTG_freqs <- c(A = 0.2487819, C = 0.2512324, T = 0.2512223, G = 0.2487634)

# Get the SOX8 motif from JASPAR2020
SOX8_tfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM", name = "SOX8"))[[1]]

# Scan for motif occurrences in the peakset
motif_hits <- matchMotifs(SOX8_tfm, NC_PPR_combined_gr, genome = BSgenome.Ggallus.UCSC.galGal6, 
                          out = "positions", p.cutoff = 1e-04, bg = ACTG_freqs)

# Subset original GRanges object to only those with a SOX8 binding motif
NC_PPR_combined_SOX8_gr <- subsetByOverlaps(NC_PPR_combined_gr, motif_hits)


################################################################################################################################
###################### CREATING GRANGES FOR plotting all Sox8+ PPR AND NC PEAKS ################################################


# Find overlapping peaks in ss4 and ss8 between NC and PPR
ss4_PPR_Sox8_Peaks <- NC_PPR_combined_SOX8_gr[NC_PPR_combined_SOX8_gr %over% ss4_PPR_all_Peaks]
ss4_NC_Sox8_Peaks <- NC_PPR_combined_SOX8_gr[NC_PPR_combined_SOX8_gr %over% ss4_NC_all_Peaks]
ss8_PPR_Sox8_Peaks <- NC_PPR_combined_SOX8_gr[NC_PPR_combined_SOX8_gr %over% ss8_PPR_all_Peaks]
ss8_NC_Sox8_Peaks <- NC_PPR_combined_SOX8_gr[NC_PPR_combined_SOX8_gr %over% ss8_NC_all_Peaks]


NC_PPR_combined_SOX8 <- unique(c(ss4_PPR_Sox8_Peaks, ss4_NC_Sox8_Peaks, ss8_PPR_Sox8_Peaks, ss8_NC_Sox8_Peaks)) #Checking that the combined peaks equal the same as the original combine GR

# Combining PPR and NC peaks from both stages
All_Stages_PPR_Sox8 <- unique(c(ss4_PPR_Sox8_Peaks,ss8_PPR_Sox8_Peaks))
All_Stages_NC_Sox8 <- unique(c(ss4_NC_Sox8_Peaks,ss8_NC_Sox8_Peaks))

# Separating into PPR only, NC only, and shared peaks
All_Stages_PPR_only_Sox8 <- setdiff(All_Stages_PPR_Sox8, All_Stages_NC_Sox8)
All_Stages_PPR_and_NC_Sox8 <- intersect(All_Stages_PPR_Sox8, All_Stages_NC_Sox8)
All_Stages_NC_only_Sox8 <- setdiff(All_Stages_NC_Sox8, All_Stages_PPR_Sox8)

# Add metadata columns to each GRanges object
mcols(All_Stages_PPR_only_Sox8)$group <- "PPR"
mcols(All_Stages_PPR_and_NC_Sox8)$group <- "Shared"
mcols(All_Stages_NC_only_Sox8)$group <- "NC"

# Now combine the ranges, placing the overlaps in the middle
Combined_SOX8_gr <- c(All_Stages_PPR_only_Sox8, All_Stages_PPR_and_NC_Sox8, All_Stages_NC_only_Sox8)
sort(Combined_SOX8_gr[ranges(Combined_SOX8_gr) %in% ranges(Combined_SOX8_gr)[duplicated(ranges(Combined_SOX8_gr))]]) # To check and return any duplicates


# Function to compute and aggregate peak accessibility scores as z-scores for heatmaps
get_enrichment_matrix_z_scores <- function(proj, peak_matrix, peaks, cell_types) {
 
  # Extract peak accessibility matrix
  mat_norm <- assay(peak_matrix)
  
   # Normalize before subsetting
  #mat_norm <- sweep(mat, 2, colSums(mat), "/") 
  
  # Retain only peaks present in "peaks" GRanges object (subset matrix first)
  overlaps <- findOverlaps(peaks, rowRanges(peak_matrix))
  matching_indices <- queryHits(overlaps)
  mat_norm <- mat_norm[matching_indices, , drop = FALSE]
  
  # Extract the matched GRanges and format as "chr:start-end"
  matched_peaks <- peaks[matching_indices]
  peak_coords <- paste0(seqnames(matched_peaks), ":", start(matched_peaks), "-", end(matched_peaks))
  
  # Extract metadata for cell types
  meta <- as.data.frame(getCellColData(proj))
  
  # Subset to relevant cell types
  selected_cells <- meta$transferred_scHelper_cell_type_broad %in% cell_types
  mat_norm <- mat_norm[, selected_cells, drop = FALSE]
  meta_filtered <- meta[selected_cells, ]
  
  # Compute row-wise z-score
  mat_scaled <- t(apply(mat_norm, 1, scale))
  rownames(mat_scaled) <- peak_coords  # Assign peak coordinates as row names
  
  # Compute mean z-score per cell type
  cell_type_labels <- meta_filtered$transferred_scHelper_cell_type_broad
  avg_mat <- t(apply(mat_scaled, 1, function(row) tapply(row, cell_type_labels, mean, na.rm = TRUE)))
  
  rownames(avg_mat) <- peak_coords  # Ensure peak coordinates remain as row names
  return(avg_mat)
}

# Function to subset and binarize the peak matrix (reads only range between 0 and 4) and hen generate a proportion of the NC and PPR cells that have at least 1 read per peak
get_enrichment_matrix_read_proporation <- function(proj, peak_matrix, peaks, cell_types) {
  
  # Extract peak accessibility matrix
  mat <- assay(peak_matrix)
  
  # Retain only peaks present in "peaks" GRanges object (subset matrix first)
  overlaps <- findOverlaps(peaks, rowRanges(peak_matrix))
  matching_indices <- queryHits(overlaps)
  mat <- mat[matching_indices, , drop = FALSE]
  
  # Extract the matched GRanges and format as "chr:start-end"
  matched_peaks <- peaks[matching_indices]
  peak_coords <- paste0(seqnames(matched_peaks), ":", start(matched_peaks), "-", end(matched_peaks))
  
  # Extract metadata for cell types
  meta <- as.data.frame(getCellColData(proj))
  
  # Subset to relevant cell types
  selected_cells <- meta$transferred_scHelper_cell_type_broad %in% cell_types
  mat <- mat[, selected_cells, drop = FALSE]
  meta_filtered <- meta[selected_cells, ]
  
  # Binarize the matrix (convert all nonzero values to 1)
  mat_binary <- mat
  mat_binary[mat > 0] <- 1
  
  # Convert binary matrix to sparse format
  mat_sparse <- as(mat_binary, "CsparseMatrix")
  
  # Create binary indicator vectors for NC and Placodal cell types (1 for TRUE, 0 for FALSE)
  NC_cells <- as.numeric(meta_filtered$transferred_scHelper_cell_type_broad == "NC")
  Placodal_cells <- as.numeric(meta_filtered$transferred_scHelper_cell_type_broad == "Placodal")
  
  # Convert the numeric vectors to sparse format using sparseVector from Matrix package
  NC_indicator <- sparseVector(x = NC_cells[NC_cells != 0], i = which(NC_cells != 0), length = length(NC_cells))
  Placodal_indicator <- sparseVector(x = Placodal_cells[Placodal_cells != 0], i = which(Placodal_cells != 0), length = length(Placodal_cells))
  
  # Calculate the sum of 1's for each row (each peak) for NC and Placodal cell types
  sum_NC <- mat_sparse %*% NC_indicator
  sum_Placodal <- mat_sparse %*% Placodal_indicator
  
  # Convert the results to a vector (remove sparse matrix structure)
  sum_NC <- as.vector(sum_NC)
  sum_Placodal <- as.vector(sum_Placodal)
  
  # Get the total number of NC and Placodal cells
  num_NC <- sum(NC_cells)
  num_Placodal <- sum(Placodal_cells)
  
  # Calculate the proportions (proportion of NC/Placodal cells with at least one read per peak)
  proportions_NC <- sum_NC / num_NC
  proportions_Placodal <- sum_Placodal / num_Placodal
  
  # Combine the proportions into a matrix
  proportions <- cbind(proportions_NC, proportions_Placodal)
  
  # Assign peak coordinates as row names
  rownames(proportions) <- peak_coords
  
  return(proportions)
}



# Define cell types of interest
cell_types <- c("Placodal", "NC")

ss4_peak_matrix <- getMatrixFromProject(ss4ArchRProj, useMatrix = "PeakMatrix")
ss8_peak_matrix <- getMatrixFromProject(ss8ArchRProj, useMatrix = "PeakMatrix")

# Get enrichment matrices for each developmental stage
mat_ss4_z_scores <- get_enrichment_matrix_z_scores(ss4ArchRProj, ss4_peak_matrix, Combined_SOX8_gr, cell_types)
mat_ss8_z_scores <- get_enrichment_matrix_z_scores(ss8ArchRProj, ss8_peak_matrix, Combined_SOX8_gr, cell_types)

mat_ss4_access_co <- get_enrichment_matrix_read_proporation(ss4ArchRProj, ss4_peak_matrix, Combined_SOX8_gr, cell_types)
mat_ss8_access_co <- get_enrichment_matrix_read_proporation(ss8ArchRProj, ss8_peak_matrix, Combined_SOX8_gr, cell_types)

# Combine both matrices for a single heatmap
mat_combined <- cbind(mat_ss4_z_scores, mat_ss8_z_scores)
range(mat_combined)

# Create column annotations (for cell type and stage)
col_anno <- data.frame(
  Stage = rep(c("ss4", "ss8"), each = length(cell_types)),
  CellType = rep(cell_types, times = 2)
)
col_anno <- HeatmapAnnotation(df = col_anno, col = list(
  Stage = c("ss4" = "lightgreen", "ss8" = "darkgreen"),
  CellType = c("Placodal" = "blue", "NC" = "red")
))

# Plot heatmap
Heatmap(mat_combined,
        name = "PPR_and_NC_peak_accessibilities",
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation = col_anno)

write.table(mat_combined, "/data/Sox8_binding_partner_analysis/scATACseq_objects/plots/Heatmaps/Sox8_Heatmap_matrix_norm_pre_subsetting")




















### troubleshooting function
  
# Extract peak accessibility matrix
# Extract peak accessibility matrix
# Extract peak accessibility matrix
mat <- assay(ss8_peak_matrix)

# Retain only peaks present in "peaks" GRanges object (subset matrix first)
overlaps <- findOverlaps(Combined_SOX8_gr, rowRanges(ss8_peak_matrix))
matching_indices <- queryHits(overlaps)
mat <- mat[matching_indices, , drop = FALSE]

# Extract the matched GRanges and format as "chr:start-end"
matched_peaks <- Combined_SOX8_gr[matching_indices]
peak_coords <- paste0(seqnames(matched_peaks), ":", start(matched_peaks), "-", end(matched_peaks))

# Extract metadata for cell types
meta <- as.data.frame(getCellColData(ss8ArchRProj))

# Subset to relevant cell types
selected_cells <- meta$transferred_scHelper_cell_type_broad %in% cell_types
mat <- mat[, selected_cells, drop = FALSE]
meta_filtered <- meta[selected_cells, ]

# Binarize the matrix (convert all nonzero values to 1)
mat_binary <- mat
mat_binary[mat > 0] <- 1

# Convert binary matrix to sparse format
mat_sparse <- as(mat_binary, "CsparseMatrix")

# Create binary indicator vectors for NC and Placodal cell types (1 for TRUE, 0 for FALSE)
NC_cells <- as.numeric(meta_filtered$transferred_scHelper_cell_type_broad == "NC")
Placodal_cells <- as.numeric(meta_filtered$transferred_scHelper_cell_type_broad == "Placodal")

# Convert the numeric vectors to sparse format using sparseVector from Matrix package
NC_indicator <- sparseVector(x = NC_cells[NC_cells != 0], i = which(NC_cells != 0) - 1, length = length(NC_cells))
Placodal_indicator <- sparseVector(x = Placodal_cells[Placodal_cells != 0], i = which(Placodal_cells != 0), length = length(Placodal_cells))

# Calculate the sum of 1's for each row (each peak) for NC and Placodal cell types
sum_NC <- mat_sparse %*% NC_indicator
sum_Placodal <- mat_sparse %*% Placodal_indicator

# Convert the results to a vector (remove sparse matrix structure)
sum_NC <- as.vector(sum_NC)
sum_Placodal <- as.vector(sum_Placodal)

# Get the total number of NC and Placodal cells
num_NC <- sum(NC_cells)
num_Placodal <- sum(Placodal_cells)

# Calculate the proportions (proportion of NC/Placodal cells with at least one read per peak)
proportions_NC <- sum_NC / num_NC
proportions_Placodal <- sum_Placodal / num_Placodal

# Combine the proportions into a matrix
proportions <- cbind(proportions_NC, proportions_Placodal)

# Assign peak coordinates as row names
rownames(proportions) <- peak_coords

return(proportions)
