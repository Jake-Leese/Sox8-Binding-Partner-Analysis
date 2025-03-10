# Script to subset an archR project peakset, in this case based on the presence of a Sox8 binding motif

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

source("/data/Sox8_binding_partner_analysis/Functions/ExtractArchRPeakMatrix.R")

## Creating genome and gene annotations for ArchR to interact with

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Ggallus.UCSC.galGal6)

geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Ggallus.UCSC.galGal6.refGene, 
                                       OrgDb = org.Gg.eg.db)

# loading ArchR projects
ss4ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Save-ArchR', force = FALSE)
ss8ArchRProj <- loadArchRProject(path = '/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Save-ArchR', force = FALSE)

# Read in existing Peak X Motif matrices
ss4_vert_motif_matrix <- readRDS("MotifMatrices/SS4_vert_motif_matrix_1e-05.rds")
ss4_Hs_motif_matrix <- readRDS("MotifMatrices/SS4_Hs_motif_matrix_1e-05.rds")
ss8_vert_motif_matrix <- readRDS("MotifMatrices/SS8_vert_motif_matrix_1e-05.rds")
ss8_Hs_motif_matrix <- readRDS("MotifMatrices/SS8_Hs_motif_matrix_1e-05.rds")


### SUBSETTING ARCHR PROJECTS BY SOX8 POSITIVE PEAKS ###

# Start by selecting only the peaks with a Sox8 binding site. Peakset is the same across both ArchR projects as it must have been generated using the full dataset.
# Sox8 binding motif matrix must be the same across vertebrate-wide and Homo sapiens-specific JASPAR2024 databases. Same number of hits. Higher hit rate than from JASPAR2020 (18266 vs 3002 peaks)
ss4_vert_Sox8_pos_peaks <- ss4_vert_motif_matrix[ss4_vert_motif_matrix[,"SOX8"]==TRUE,c(1:879)]
dim(ss4_vert_Sox8_pos_peaks) # 18266 X 879
ss4_Hs_Sox8_pos_peaks <- ss4_Hs_motif_matrix[ss4_Hs_motif_matrix[,"SOX8"]==TRUE,c(1:720)]
dim(ss4_Hs_Sox8_pos_peaks)   # 18266 X 879


# Extract the seqnames (chr) and ranges (Iranges) from the Sox8_pos_peaks 

Sox8_peaks <- (rownames(ss4_Hs_Sox8_pos_peaks)) # Extract Genome co-ordinates of Sox8 peaks
DF <- data.frame(Sox8_peaks)                    # Create a dataframe, DF, that separates seqnames (chrx) from genomic ranges.
DF <- DF %>%                                    # Issue with this, is that it separates out values by the presence of a "-"
  separate(                                     # Ranges column therefore only contains the first value of the IRange before the "-"
    Sox8_peaks,
    into = c("seqnames", "ranges"),
    sep = "-",
    extra = "drop")
Sox8_seqnames <- as.vector(DF$seqnames)

# Function to extract characters after the first dash only
extract_after_first_dash <- function(input_vector) {
  # Initialize an empty vector to store the results
  result_vector <- character(length(input_vector))
  
  # Loop through each element in the input vector
  for (i in seq_along(input_vector)) {
    # Find the position of the first "-"
    first_dash_position <- regexpr("-", input_vector[i])
    
    # Check if a dash is found
    if (first_dash_position != -1) {
      # Extract the substring starting from the position after the first "-"
      result_vector[i] <- substr(input_vector[i], first_dash_position + 1, nchar(input_vector[i]))
    } else {
      # If no dash is found, retain the entire string
      result_vector[i] <- input_vector[i]
    }
  }
  
  return(result_vector)
}

Sox8_ranges <- extract_after_first_dash(Sox8_peaks)

# Create a new GRanges object using the seqnames and ranges
Sox8_GR <- GRanges(seqnames = Sox8_seqnames, ranges = Sox8_ranges)


#### SUBSETTING ARCHR PROJECT BASED ON NEWLY MADE SOX8 PEAK GRANGES OBJECT ####

# First we need to extract the peakset from ArchR project. (can only have 1 peakset per archrproj)
ss8_Full_peakSet <- getPeakSet(ArchRProj = ss8ArchRProj)

# Subset the ss8 ArchR project Peakset using subsetByOverlaps() function. Save to disk for future use.
# This is almost identical to the Sox8_GR we just made, except that it retains the metadata from the ArchR project, which we want to keep when we come to add it back in
ss8_Sox8_positive_peakset <- subsetByOverlaps(ss8_Full_peakSet, Sox8_GR)
saveRDS(ss8_Sox8_positive_peakset, file = "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/SS8_Sox8_positive_peakset.rds")

# Adding new subsetted peakset to a new ArchR project, and saving this new ArchR project.
ss8_Sox8_ArchRProj <- addPeakSet(
  ArchRProj = ss8ArchRProj,
  peakSet = ss8_Sox8_positive_peakset,
  genomeAnnotation = getGenomeAnnotation(ss8ArchRProj),
  force = TRUE)

# Adding PeakMatrix using new Sox8_positive_peakset
ss8_Sox8_ArchRProj <- addPeakMatrix(ss8_Sox8_ArchRProj)

#Save ArchRProject
ss8_Sox8_ArchRProj <- saveArchRProject(ss8_Sox8_ArchRProj, outputDirectory = "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss8_Sox8_Save_ArchR", load = TRUE)
ss8_Sox8_ArchRProj@peakSet  # Check peakset

### Reapeating all of above with ss4 ArchR proj
ss4_Full_peakSet <- getPeakSet(ArchRProj = ss4ArchRProj)

ss4_Sox8_positive_peakset <- subsetByOverlaps(ss4_Full_peakSet, Sox8_GR)
saveRDS(ss4_Sox8_positive_peakset, file = "/data/Sox8_binding_partner_analysis/scATACseq_objects/Peaksets/SS4_Sox8_positive_peakset.rds")

ss4_Sox8_ArchRProj <- addPeakSet(
  ArchRProj = ss4ArchRProj,
  peakSet = ss4_Sox8_positive_peakset,
  genomeAnnotation = getGenomeAnnotation(ss4ArchRProj),
  force = TRUE)

ss4_Sox8_ArchRProj <- addPeakMatrix(ss4_Sox8_ArchRProj)

ss4_Sox8_ArchRProj <- saveArchRProject(ss4_Sox8_ArchRProj, outputDirectory = "/data/Sox8_binding_partner_analysis/scATACseq_objects/ss4_Sox8_Save_ArchR", load = TRUE)
ss4_Sox8_ArchRProj@peakSet


#### COMPARING THE PEAK MATRICES BETWEEN ORIGINAL AND SUBSETTED ARCHR PROJECTS ####

ss4_full_peakMatrix <- ExtractArchRPeakMatrix(ss4ArchRProj)       # Extract labelled peak matrices using custom function
ss4_Sox8_peakMatrix <- ExtractArchRPeakMatrix(ss4_Sox8_ArchRProj)

common_rows <- intersect(rownames(ss4_full_peakMatrix), rownames(ss4_Sox8_peakMatrix))    # Extract rownames common to both peak matrices
common_columns <- intersect(colnames(ss4_full_peakMatrix), colnames(ss4_Sox8_peakMatrix)) # Extract colnames common to both peak matrices

sub_ss4_full_peakMatrix <- ss4_full_peakMatrix[common_rows, common_columns] # Subset peak matrices based on common row and col names
sub_ss4_Sox8_peakMatrix <- ss4_Sox8_peakMatrix[common_rows, common_columns]

comparison <- sub_ss4_full_peakMatrix - sub_ss4_Sox8_peakMatrix # Compare values of subsetted matrices, which should now have the same peak and cell labels
comparison[comparison > 0] # =0 therefore the peak matrix values are the same between both archR projects
