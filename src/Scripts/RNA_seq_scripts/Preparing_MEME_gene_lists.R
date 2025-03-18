# Script to follow the same process used to characterise Sox8 co-expressed TFs, but extending it to all genes
# These will then be used to make gene lists for enhancer annotation of Placodal and NC enhancers
# Run using ArchR_Seurat_R_4.4.1.sif container

setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(TFBSTools)
library(JASPAR2024)
library(dplyr)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(universalmotif)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)

############################### MAKE PPR/NC GENE LISTS FROM SCRNA-SEQ DATA ##########################################

source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_HM_Clusters.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/subset_seurat_features.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Extract_count_data.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/Correlation_analysis.R")
source("/data/Sox8_binding_partner_analysis/src/Functions/plot_correlation_against_counts.R")


# Load the original Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")


####################################################################################################################
######################## DEFINING GENE EXPRESSION BASED ON CORRELATION WITH SOX8 ###################################


# Custom function to extract count data for: Stage, Cell_type (ectoderm type), Sample (origin.ident)(optional), and Assay ("RNA" or "SCT"), data_type ("counts" or "data" (log.normalised)), and subset_cells_by_gene ("gene.name") 
placode_HH8_RNA_Counts <- Extract_count_data(HHall_ectoderm, Stage = "HH8", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
NC_HH8_RNA_Counts <- Extract_count_data(HHall_ectoderm, Stage = "HH8", Cell_type = "NC", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
placode_HH9_RNA_Counts <- Extract_count_data(HHall_ectoderm, Stage = "HH9", Cell_type = "placode", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")
NC_HH9_RNA_Counts <- Extract_count_data(HHall_ectoderm, Stage = "HH9", Cell_type = "NC", Assay = "RNA", data_type = "counts", subset_cells_by_gene = "SOX8")

#### CALCULATE CORRELATION COEFFICIENTS ####
# Returns correlation values for specified gene only
placode_HH8_RNA_spearman_SOX8 <- Correlation_analysis(placode_HH8_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
NC_HH8_RNA_spearman_SOX8 <- Correlation_analysis(NC_HH8_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
placode_HH9_RNA_spearman_SOX8 <- Correlation_analysis(placode_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")
NC_HH9_RNA_spearman_SOX8 <- Correlation_analysis(NC_HH9_RNA_Counts, remove_null_count_genes = TRUE, Test = "spearman", Gene = "SOX8")


# extract values from gene specific correlation dataframe and format appropriately for plotting
placode_HH8_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(placode_HH8_RNA_spearman_SOX8[placode_HH8_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH8_placode") 
NC_HH8_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(NC_HH8_RNA_spearman_SOX8[NC_HH8_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH8_NC")
placode_HH9_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(placode_HH9_RNA_spearman_SOX8[placode_HH9_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH9_placode")
NC_HH9_RNA_spearman_cor_val <- data.frame(correlation = as.numeric(NC_HH9_RNA_spearman_SOX8[NC_HH9_RNA_spearman_SOX8 < 0.99]), Stage_CellType = "HH9_NC")

# combine dataframes 
Combined_df <- rbind(placode_HH8_RNA_spearman_cor_val, NC_HH8_RNA_spearman_cor_val, placode_HH9_RNA_spearman_cor_val, NC_HH9_RNA_spearman_cor_val)

# Calculate the mean + SD for each group
#thresholds <- Combined_df %>%
 #group_by(Stage_CellType) %>%
  #summarize(threshold = mean(correlation) + sd(correlation))

# Filter genes based on thresholds
#HH8_Placode_genes <- colnames(placode_HH8_RNA_spearman_SOX8)[placode_HH8_RNA_spearman_SOX8 >= 0.103]   # 17426 to 2773
#HH8_NC_genes <- colnames(NC_HH8_RNA_spearman_SOX8)[NC_HH8_RNA_spearman_SOX8 >= 0.0994]                 # 17452 to 2594
#HH9_Placode_genes <- colnames(placode_HH9_RNA_spearman_SOX8)[placode_HH9_RNA_spearman_SOX8 >= 0.159]   # 17568 to 2905
#HH9_NC_genes <- colnames(NC_HH9_RNA_spearman_SOX8)[NC_HH9_RNA_spearman_SOX8 >= 0.233]                  # 21379 to 3477

# Filter based on previously calculated thresholds (those used for transcription factors)
HH8_Placode_genes <- colnames(placode_HH8_RNA_spearman_SOX8)[placode_HH8_RNA_spearman_SOX8 >= 0.109216]   # 17426 to 2773
HH8_NC_genes <- colnames(NC_HH8_RNA_spearman_SOX8)[NC_HH8_RNA_spearman_SOX8 >= 0.117146]                 # 17452 to 2594
HH9_Placode_genes <- colnames(placode_HH9_RNA_spearman_SOX8)[placode_HH9_RNA_spearman_SOX8 >= 0.148648]   # 17568 to 2905
HH9_NC_genes <- colnames(NC_HH9_RNA_spearman_SOX8)[NC_HH9_RNA_spearman_SOX8 >= 0.208392]                  # 21379 to 3477

############################### ASSIGN GALGAL6 ACCESSION NUMBERS TO TFS 
galgal6_gtf <- import("/data/Sox8_binding_partner_analysis/genome_files/Gallus_gallus.GRCg6a.97.gtf")
galgal6_gtf_df <- as.data.frame(galgal6_gtf)
  
# Function to return character vector of gene names, gene_ids, or gene_list, and return gene name and corresponding gene_ids
Add_gene_IDs <- function(gtf, gene_list) {
filtered_gtf <- gtf %>%
  filter(str_detect(tolower(gene_name), tolower(paste(gene_list, collapse = "|"))) | gene_id %in% gene_list | transcript_id %in% gene_list) %>%
  select(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)
}  

HH8_Placode_gene_ids <- Add_gene_IDs(galgal6_gtf_df, HH8_Placode_genes)
HH8_NC_gene_ids <- Add_gene_IDs(galgal6_gtf_df, HH8_NC_genes)
HH9_Placode_gene_ids <- Add_gene_IDs(galgal6_gtf_df, HH9_Placode_genes)
HH9_NC_gene_ids <- Add_gene_IDs(galgal6_gtf_df, HH9_NC_genes)

# Adding Sox8 coexpression values to each gene for future reference
HH8_Placode_gene_ids$Sox8_coexpression <- t(placode_HH8_RNA_spearman_SOX8[HH8_Placode_gene_ids$gene_name])
HH8_NC_gene_ids$Sox8_coexpression <- t(NC_HH8_RNA_spearman_SOX8[HH8_NC_gene_ids$gene_name])
HH9_Placode_gene_ids$Sox8_coexpression <- t(placode_HH9_RNA_spearman_SOX8[HH9_Placode_gene_ids$gene_name])
HH9_NC_gene_ids$Sox8_coexpression <- t(NC_HH9_RNA_spearman_SOX8[HH9_NC_gene_ids$gene_name])

# Saving as .csv file
write.csv(HH8_Placode_gene_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_Placode_gene_ids.csv")
write.csv(HH8_NC_gene_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_gene_ids.csv")
write.csv(HH9_Placode_gene_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_Placode_gene_ids.csv")
write.csv(HH9_NC_gene_ids, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_gene_ids.csv")

# Saving as gene_lists using gene_ids for enhancer annotation
writeLines(HH8_Placode_gene_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_Placode_genes.txt")
writeLines(HH8_NC_gene_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_genes.txt")
writeLines(HH9_Placode_gene_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_Placode_genes.txt")
writeLines(HH9_NC_gene_ids$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_NC_genes.txt")

# Saving unfiltered gene_lists using default names from seurat
#writeLines(HH8_Placode_genes, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_Placode_genes.txt")
#writeLines(HH8_NC_genes, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_genes.txt")
#writeLines(HH9_Placode_genes, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_Placode_genes.txt")
#writeLines(HH9_NC_genes, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_NC_genes.txt")


####################################################################################################################
########################## DEFINING SEURAT DEGS AS DOWNSTREAM GENES OF INTEREST ####################################

# Subsetting seurat object cells by stage
HH8_ectoderm_genes <- HHall_ectoderm[,HHall_ectoderm$stage %in% "HH8"]
HH9_ectoderm_genes <- HHall_ectoderm[,HHall_ectoderm$stage %in% "HH9"]

# Perform differential expression analysis between NC and PPR for each stage
HH8_Placodal_vs_NC <- FindMarkers(
  HH8_ectoderm_genes,
  group.by = "ectoderm_type",
  ident.1 = "placode",
  ident.2 = "NC",
  assay = "SCT",
  min.pct = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)

HH9_Placodal_vs_NC <- FindMarkers(
  HH9_ectoderm_genes,
  group.by = "ectoderm_type",
  ident.1 = "placode",
  ident.2 = "NC",
  assay = "SCT",
  min.pct = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)

# Extract genes expressed significantly more in Placode or NC cells. + log2FC is placodal, - log2FC is NC
HH8_PPR_DEGs <- rownames(HH8_Placodal_vs_NC[HH8_Placodal_vs_NC$avg_log2FC > 0 & HH8_Placodal_vs_NC$p_val_adj < 0.05, ])
HH8_NC_DEGs <- rownames(HH8_Placodal_vs_NC[HH8_Placodal_vs_NC$avg_log2FC < 0 & HH8_Placodal_vs_NC$p_val_adj < 0.05, ])
HH9_PPR_DEGs <- rownames(HH9_Placodal_vs_NC[HH9_Placodal_vs_NC$avg_log2FC > 0 & HH9_Placodal_vs_NC$p_val_adj < 0.05, ])
HH9_NC_DEGs <- rownames(HH9_Placodal_vs_NC[HH9_Placodal_vs_NC$avg_log2FC < 0 & HH9_Placodal_vs_NC$p_val_adj < 0.05, ])


# Adding gene_ids from GalGal6 GTF using function defined in previous section
HH8_PPR_DEGs_IDs <- Add_gene_IDs(galgal6_gtf_df, HH8_PPR_DEGs)
HH8_NC_DEGs_IDs <- Add_gene_IDs(galgal6_gtf_df, HH8_NC_DEGs)
HH9_PPR_DEGs_IDs <- Add_gene_IDs(galgal6_gtf_df, HH9_PPR_DEGs)
HH9_NC_DEGs_IDs <- Add_gene_IDs(galgal6_gtf_df, HH9_NC_DEGs)

# Saving as .csv file
write.csv(HH8_PPR_DEGs_IDs, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_Placode_DEG_ids_gal6.csv")
write.csv(HH8_NC_DEGs_IDs, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_DEG_ids_gal6.csv")
write.csv(HH9_PPR_DEGs_IDs, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_Placode_DEG_ids_gal6.csv")
write.csv(HH9_NC_DEGs_IDs, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_DEG_ids_gal6.csv")

# Saving as gene_lists using gene_ids for enhancer annotation
writeLines(HH8_PPR_DEGs_IDs$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_Placode_DEGs_gal6.txt")
writeLines(HH8_NC_DEGs_IDs$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_DEGs_gal6.txt")
writeLines(HH9_PPR_DEGs_IDs$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_Placode_DEGs_gal6.txt")
writeLines(HH9_NC_DEGs_IDs$gene_id, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_NC_DEGs_gal6.txt")



##########################################################################################################################
##################################### MATCHING GALGAL7 TO GALGAL6 GENES ##################################################

# Comparing GalGal6 and GalGal7 gene names
galgal6_all_gene_names <- unique(galgal6_gtf_df$gene_name) # 13,486 unique gene names
galgal7_all_gene_names <- unique(galgal7_gtf_df$gene_name) # 13,975 unique gene names
length(intersect(galgal6_all_gene_names, galgal7_all_gene_names)) #12,093 common gene names across both (out of 13,486 galgal6 genes)
length(intersect(tolower(galgal6_all_gene_names), tolower(galgal7_all_gene_names))) # 12,106 when we remove case sensitivity
sum(str_detect(tolower(galgal6_all_gene_names), tolower(paste(galgal7_all_gene_names, collapse = "|"))), na.rm = TRUE) # 12,539 when we count all cases of galgal7 gene names whose string can be detected within galgal6 gene names, even if they don't perfectly match

# Checking intersect between gene_names in our galgal7 DEG lists and those that are not found (setdiff) in GalGal6
HH8_PPR_DEGs_GalGal7_only <- intersect(HH8_PPR_DEGs_IDs$gene_name, setdiff(galgal7_all_gene_names, galgal6_all_gene_names))  # 105 "unique" galgal7 genes
HH8_NC_DEGs_GalGal7_only <- intersect(HH8_NC_DEGs_IDs$gene_name, setdiff(galgal7_all_gene_names, galgal6_all_gene_names))  # 83 "unique" galgal7 genes
HH9_PPR_DEGs_GalGal7_only <- intersect(HH9_PPR_DEGs_IDs$gene_name, setdiff(galgal7_all_gene_names, galgal6_all_gene_names))  # 188 "unique" galgal7 genes
HH9_NC_DEGs_GalGal7_only <- intersect(HH9_NC_DEGs_IDs$gene_name, setdiff(galgal7_all_gene_names, galgal6_all_gene_names))  # 202 "unique" galgal7 genes


# GalGal6 genes that are saved by removing case sensitivity and using str_detect() in Add_gene_IDs()
HH8_PPR_galgal6_recovered <- HH8_PPR_DEGs_IDs[str_detect(tolower(HH8_PPR_DEGs_IDs$gene_name), tolower(paste(HH8_PPR_DEGs_GalGal7_only, collapse = "|"))), ] # Recovers 33 of the 105 "unique" galgal7 genes
HH8_NC_galgal6_recovered <- HH8_NC_DEGs_IDs[str_detect(tolower(HH8_NC_DEGs_IDs$gene_name), tolower(paste(HH8_NC_DEGs_GalGal7_only, collapse = "|"))), ] # Recovers 16 of the 83 "unique" galgal7 genes
HH9_PPR_galgal6_recovered <- HH9_PPR_DEGs_IDs[str_detect(tolower(HH9_PPR_DEGs_IDs$gene_name), tolower(paste(HH9_PPR_DEGs_GalGal7_only, collapse = "|"))), ] # Recovers 43 of the 188 "unique" galgal7 genes
HH9_NC_galgal6_recovered <- HH9_NC_DEGs_IDs[str_detect(tolower(HH9_NC_DEGs_IDs$gene_name), tolower(paste(HH9_NC_DEGs_GalGal7_only, collapse = "|"))), ]  # Recovers 26 of the 202 "unique" galgal7 genes

# Saving matched genes to file for future reference
writeLines(HH8_PPR_galgal6_recovered$gene_name , "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_PPR_DEGs_gal6_recovered.txt")
writeLines(HH8_NC_galgal6_recovered$gene_name , "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH8_NC_DEGs_gal6_recovered.txt")
writeLines(HH9_PPR_galgal6_recovered$gene_name , "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_PPR_DEGs_gal6_recovered.txt")
writeLines(HH9_NC_galgal6_recovered$gene_name , "/data/Sox8_binding_partner_analysis/scRNAseq_objects/CoExpressed_genes/HH9_NC_DEGs_gal6_recovered.txt")


# Searching gtf dataframes for specific gene names and IDs
galgal6_gtf_df[galgal6_gtf_df$gene_id == "ENSGALG00000039118", ]
galgal7_gtf_df_no_na <- galgal7_gtf_df[!is.na(galgal7_gtf_df$gene_name), ]
galgal7_gtf_df_no_na[galgal7_gtf_df_no_na$gene_name == "MEIS2", ]

