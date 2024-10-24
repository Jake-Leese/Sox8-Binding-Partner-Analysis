setwd("/data/Sox8_binding_partner_analysis/scRNAseq_objects/")
.libPaths("/R/libs/ArchR_Seurat_R_441")

library(Seurat)
library(usethis)
library(devtools)
library(CSCORE)
library(TFBSTools)
library(JASPAR2024)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(universalmotif)

source("/data/Sox8_binding_partner_analysis/Functions/Extract_HM_Clusters.R")

# Load and check the Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

head(HHall_ectoderm)


### SUBSET COUNT DATA FOR JASPAR2024 TFs ONLY ####

# Download JASPAR2024 motif list for vertebrates
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

# extract list of genes from Seurat object and checking overlap with TF lists
# RNAseq_all_genes <- rownames(HHall_ectoderm)
# vert_tf_genes <- intersect(RNAseq_all_genes, name_vector_vert) # 422 overlapping genes
# Hs_tf_genes <- intersect(RNAseq_all_genes, name_vector_Hs)     # 417 overlapping genes 
# setdiff(vert_tf_genes, Hs_tf_genes) # Genes missing from Hs_ts_genes: TBP, LIN54, NR5A1, Dlx4, Atoh1. Will move forward with vertebrate TFs


# Subsetting Seurat object to include JASPAR TFs only
HHall_ectoderm_TFs <- DietSeurat(HHall_ectoderm, features = name_vector_vert)


#### EXTRACT COUNT DATA FOR STAGES HH8-HH9 AND CELL TYPES NC AND PLACODAL ####

# Extract and subset cells based on stage (HH8-HH9) and cell type (NC and Placodal)
HHall_NC_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$ectoderm_type %in% "NC"]
HH8_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH8"]
HH9_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH9"]

HHall_Placodal_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$ectoderm_type %in% "placode"]
HH8_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH8"]
HH9_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH9"]

# extract count data from NC_only and placode_only SEURAT object, and then subset based on TF genes only
HH8_NC_TFs_counts <- t(HH8_NC_TFs[["SCT"]]$counts)
HH8_NC_TFs_counts
HH9_NC_TFs_counts <- t(HH9_NC_TFs[["SCT"]]$counts)
HH8_Placodal_TFs_counts <- t(HH8_Placodal_TFs[["SCT"]]$counts)
HH9_Placodal_TFs_counts <- t(HH9_Placodal_TFs[["SCT"]]$counts)

# calculate correlation matrices using spearman
HH8_NC_TF_correlation_matrix <- WGCNA::cor(HH8_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH8_NC <- data.frame(as.list(HH8_NC_TF_correlation_matrix["SOX8", ]))

HH9_NC_TF_correlation_matrix <- WGCNA::cor(HH9_NC_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH9_NC <- data.frame(as.list(HH9_NC_TF_correlation_matrix["SOX8", ]))

HH8_Placodal_TF_correlation_matrix <- WGCNA::cor(HH8_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH8_Placodal <- data.frame(as.list(HH8_Placodal_TF_correlation_matrix["SOX8", ]))

HH9_Placodal_TF_correlation_matrix <- WGCNA::cor(HH9_Placodal_TFs_counts, method = "spearman")
Sox8_co_expressed_TFs_HH9_Placodal <- data.frame(as.list(HH9_Placodal_TF_correlation_matrix["SOX8", ]))





# Filter out genes that have a co-expression coefficient of 0.1 or greater
HH8_NC_Sox8_coexpression_0.1 <- Sox8_co_expressed_TFs_HH8_NC[which(Sox8_co_expressed_TFs_HH8_NC > 0.1)]
HH8_Placodal_Sox8_coexpression_0.1 <- Sox8_co_expressed_TFs_HH8_Placodal[which(Sox8_co_expressed_TFs_HH8_Placodal > 0.1)]
HH9_NC_Sox8_coexpression_0.1 <- Sox8_co_expressed_TFs_HH9_NC[which(Sox8_co_expressed_TFs_HH9_NC > 0.1)]
HH9_Placodal_Sox8_coexpression_0.1 <- Sox8_co_expressed_TFs_HH9_Placodal[which(Sox8_co_expressed_TFs_HH9_Placodal > 0.1)]

HH8_NC_only <- setdiff(HH8_NC_Sox8_coexpression_0.1, HH8_Placodal_Sox8_coexpression_0.1)
HH8_Placodal_only <- setdiff(t(HH8_Placodal_Sox8_coexpression_0.1), t(HH8_NC_Sox8_coexpression_0.1))
HH8_NC_and_Placodal <- (intersect(colnames(HH8_NC_Sox8_coexpression_0.1), colnames(HH8_Placodal_Sox8_coexpression_0.1)))
HH8_NC_and_Placodal <- t(data.frame(
  HH8_NC = sapply(HH8_NC_and_Placodal, function(col) HH8_NC_Sox8_coexpression_0.1[[col]]),
  HH8_Placodal = sapply(HH8_NC_and_Placodal, function(col) HH8_Placodal_Sox8_coexpression_0.1[[col]])
))

# Make a dataframe from Sox8 Spearman's rank values for each stage
# First, assign rownames to each Sox8 correlation dataframe.
rownames(Sox8_co_expressed_TFs_HH8_NC) <- "HH8"
rownames(Sox8_co_expressed_TFs_HH9_NC) <- "HH9"

rownames(Sox8_co_expressed_TFs_HH8_Placodal) <- "HH8"
rownames(Sox8_co_expressed_TFs_HH9_Placodal) <- "HH9"

# Then combine the dataframes into 1 using rbind()
Sox8_NC_co_expressed_TFs_HH8_HH9 <- rbind(Sox8_co_expressed_TFs_HH8_NC,
                                          Sox8_co_expressed_TFs_HH9_NC)
Sox8_NC_co_expressed_TFs_HH8_HH9 <- t(Sox8_NC_co_expressed_TFs_HH8_HH9)

Sox8_Placodal_co_expressed_TFs_HH8_HH9 <- rbind(Sox8_co_expressed_TFs_HH8_Placodal,
                                                Sox8_co_expressed_TFs_HH9_Placodal)
Sox8_Placodal_co_expressed_TFs_HH8_HH9 <- t(Sox8_Placodal_co_expressed_TFs_HH8_HH9)

# Remove SOX8 from dataframe and Replace all NA values with a 0 and save to disk
Sox8_NC_co_expressed_TFs_HH8_HH9 <- subset(Sox8_NC_co_expressed_TFs_HH8_HH9, rownames(Sox8_NC_co_expressed_TFs_HH8_HH9) != "SOX8")
Sox8_NC_co_expressed_TFs_HH8_HH9[is.na(Sox8_NC_co_expressed_TFs_HH8_HH9)] <- 0

write.csv(Sox8_NC_co_expressed_TFs_HH8_HH9, "CoExpression_values/Sox8_NC_co_expressed_JASPAR_TFs_HH8_HH9.csv")
sum(Sox8_NC_co_expressed_TFs_HH8_HH9[,2] > 0.05)

Sox8_Placodal_co_expressed_TFs_HH8_HH9 <- subset(Sox8_Placodal_co_expressed_TFs_HH8_HH9, rownames(Sox8_Placodal_co_expressed_TFs_HH8_HH9) != "SOX8")
Sox8_Placodal_co_expressed_TFs_HH8_HH9[is.na(Sox8_Placodal_co_expressed_TFs_HH8_HH9)] <- 0

write.csv(Sox8_Placodal_co_expressed_TFs_HH8_HH9, "CoExpression_values/Sox8_Placodal_co_expressed_JASPAR_TFs_HH8_HH9.csv")
sum(Sox8_Placodal_co_expressed_TFs_HH8_HH9[,2] > 0.05)


# Alternatively, calculate correlation coefficients using pearsons and repeat steps above
HH8_NC_TF_correlation_matrix_pearson <- WGCNA::cor(HH8_NC_TFs_counts, method = "pearson")
HH8_NC_TF_correlation_matrix_pearson[is.na(HH8_NC_TF_correlation_matrix_pearson)] <- 0
HH8_NC_TF_correlation_matrix_pearson <- HH8_NC_TF_correlation_matrix_pearson[as.vector(colSums(HH8_NC_TF_correlation_matrix_pearson) != 0), as.vector(colSums(HH8_NC_TF_correlation_matrix_pearson) != 0)]
mean(abs(HH8_NC_TF_correlation_matrix_pearson)) # Mean correlation value of 0.0448758 (based on absolute values and with all 0 correlation values removed)
median(abs(HH8_NC_TF_correlation_matrix_pearson)) # median abs value is 0.03101766
max(HH8_NC_TF_correlation_matrix_pearson[HH8_NC_TF_correlation_matrix_pearson < 0.99]) # Maximum correlation value in whole matrix is only 0.7058475
Sox8_co_expressed_TFs_HH8_NC_pearson <- data.frame(as.list(HH8_NC_TF_correlation_matrix_pearson["SOX8", ]))

HH9_NC_TF_correlation_matrix_pearson <- WGCNA::cor(HH9_NC_TFs_counts, method = "pearson")
Sox8_co_expressed_TFs_HH9_NC_pearson <- data.frame(as.list(HH9_NC_TF_correlation_matrix_pearson["SOX8", ]))

HH8_Placodal_TF_correlation_matrix_pearson <- WGCNA::cor(HH8_Placodal_TFs_counts, method = "pearson")
Sox8_co_expressed_TFs_HH8_Placodal_pearson <- data.frame(as.list(HH8_Placodal_TF_correlation_matrix_pearson["SOX8", ]))

HH9_Placodal_TF_correlation_matrix_pearson <- WGCNA::cor(HH9_Placodal_TFs_counts, method = "pearson")
Sox8_co_expressed_TFs_HH9_Placodal_pearson <- data.frame(as.list(HH9_Placodal_TF_correlation_matrix_pearson["SOX8", ]))

# Make a dataframe from Sox8 Pearson's rank values for each stage
# First, assign rownames to each Sox8 correlation dataframe.
rownames(Sox8_co_expressed_TFs_HH8_NC_pearson) <- "HH8"
rownames(Sox8_co_expressed_TFs_HH9_NC_pearson) <- "HH9"

rownames(Sox8_co_expressed_TFs_HH8_Placodal_pearson) <- "HH8"
rownames(Sox8_co_expressed_TFs_HH9_Placodal_pearson) <- "HH9"

# Then combine the dataframes into 1 using rbind()
Sox8_NC_co_expressed_TFs_HH8_HH9_pearson <- rbind(Sox8_co_expressed_TFs_HH8_NC_pearson,
                                          Sox8_co_expressed_TFs_HH9_NC_pearson)
Sox8_NC_co_expressed_TFs_HH8_HH9_pearson <- t(Sox8_NC_co_expressed_TFs_HH8_HH9_pearson)

Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson <- rbind(Sox8_co_expressed_TFs_HH8_Placodal_pearson,
                                                Sox8_co_expressed_TFs_HH9_Placodal_pearson)
Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson <- t(Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson)

# Remove SOX8 from dataframe and Replace all NA values with a 0 and save to disk
Sox8_NC_co_expressed_TFs_HH8_HH9_pearson <- subset(Sox8_NC_co_expressed_TFs_HH8_HH9_pearson, rownames(Sox8_NC_co_expressed_TFs_HH8_HH9_pearson) != "SOX8")
Sox8_NC_co_expressed_TFs_HH8_HH9_pearson[is.na(Sox8_NC_co_expressed_TFs_HH8_HH9_pearson)] <- 0

write.csv(Sox8_NC_co_expressed_TFs_HH8_HH9_pearson, "CoExpression_values/Sox8_NC_co_expressed_JASPAR_TFs_HH8_HH9_pearson.csv")
sum(Sox8_NC_co_expressed_TFs_HH8_HH9_pearson[,2] > 0.05)

Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson <- subset(Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson, rownames(Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson) != "SOX8")
Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson[is.na(Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson)] <- 0

write.csv(Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson, "CoExpression_values/Sox8_Placodal_co_expressed_JASPAR_TFs_HH8_HH9_pearson.csv")
sum(Sox8_Placodal_co_expressed_TFs_HH8_HH9_pearson[,2] > 0.05)


#### READING IN SOX8 CO-EXPRESSION VALUES AND FORMATTING FOR CLUSTERING / HEATMAP PLOTTING

# Read in .csv files. Read column 1 (genes) as row names, rather than a data column
Sox8_NC_TF_HH8_HH9_co_expression <- read.csv("CoExpression_values/Sox8_NC_co_expressed_JASPAR_TFs_HH8_HH9.csv", row.names = 1)
Sox8_Placodal_TF_HH8_HH9_co_expression <- read.csv("CoExpression_values/Sox8_Placodal_co_expressed_JASPAR_TFs_HH8_HH9.csv", row.names = 1)

# Rename columns in dataframes to represent stage_celltype
colnames(Sox8_NC_TF_HH8_HH9_co_expression) <- paste(colnames(Sox8_NC_TF_HH8_HH9_co_expression), "NC", sep = "_")
colnames(Sox8_Placodal_TF_HH8_HH9_co_expression) <- paste(colnames(Sox8_Placodal_TF_HH8_HH9_co_expression), "Placodal", sep = "_")

# Combine dataframes and reorder columns 
correlation_matrix_HH8_HH9 <- cbind(Sox8_NC_TF_HH8_HH9_co_expression, Sox8_Placodal_TF_HH8_HH9_co_expression)
correlation_matrix_HH8_HH9 <- correlation_matrix_HH8_HH9[, c(1,3,2,4)] # re-ordering columns

# Filtering out all of the rows (genes) in the matrix that do not have at least one Sox8 Co-expression value of 0.1 or greater
correlation_matrix_HH8_HH9_0.1 <- correlation_matrix_HH8_HH9[apply(correlation_matrix_HH8_HH9, 1, function(row) any(row >= 0.1)),]


# Repeat above steps for pearson correlations
# Read in .csv files. Read column 1 (genes) as row names, rather than a data column
Sox8_NC_TF_HH8_HH9_co_expression_pearson <- read.csv("CoExpression_values/Sox8_NC_co_expressed_JASPAR_TFs_HH8_HH9_pearson.csv", row.names = 1)
Sox8_Placodal_TF_HH8_HH9_co_expression_pearson <- read.csv("CoExpression_values/Sox8_Placodal_co_expressed_JASPAR_TFs_HH8_HH9_pearson.csv", row.names = 1)

# Rename columns in dataframes to represent stage_celltype
colnames(Sox8_NC_TF_HH8_HH9_co_expression_pearson) <- paste(colnames(Sox8_NC_TF_HH8_HH9_co_expression_pearson), "NC", sep = "_")
colnames(Sox8_Placodal_TF_HH8_HH9_co_expression_pearson) <- paste(colnames(Sox8_Placodal_TF_HH8_HH9_co_expression_pearson), "Placodal", sep = "_")

# Combine dataframes and reorder columns 
correlation_matrix_HH8_HH9_pearson <- cbind(Sox8_NC_TF_HH8_HH9_co_expression_pearson, Sox8_Placodal_TF_HH8_HH9_co_expression_pearson)
correlation_matrix_HH8_HH9_pearson <- correlation_matrix_HH8_HH9_pearson[, c(1,3,2,4)] # re-ordering columns

# Filtering out all of the rows (genes) in the matrix that do not have at least one Sox8 Co-expression value of 0.1 or greater
correlation_matrix_HH8_HH9_0.1_pearson <- correlation_matrix_HH8_HH9_pearson[apply(correlation_matrix_HH8_HH9_pearson, 1, function(row) any(row >= 0.1)),]



#### PLOTTING DATA USING COMPLEX HEATMAP

# hierarchical clustering using ward.D2 method. Distances calculated using euclidean method.
HCluster <- hclust(dist(correlation_matrix_HH8_HH9_0.1_pearson, method = "euclidean"), method = "ward.D2")

# Converting dataframe to numerical matrix in order to be read by complex heatmap         
correlation_matrix_filtered_mat <- data.matrix(correlation_matrix_HH8_HH9_0.1_pearson)

# defining a column split vector. unique(col_split) then respects the ordering of the input vector (https://stackoverflow.com/questions/55799509/avoid-re-ordering-of-rows-and-columns-in-a-heatmap-r)
col_split <- rep(c("HH8", "HH9"), each = 2)

# Creating heatmap column labels and defining colours
stages <- c("HH8", "HH8", "HH9", "HH9")
CellType <- c("NC", "Placodal", "NC", "Placodal")

annotation_col <- data.frame(
  Stage = stages,
  CellType = CellType)
rownames(annotation_col) <- colnames(correlation_matrix_HH8_HH9_0.1_pearson)

# Define colors for annotation
annotation_colors <- list(
  Stage = c(HH8 = "lightblue", HH9 = "darkolivegreen2"),
  CellType = c(NC = "sienna3", Placodal = "darkcyan")
)

# Plotting complex heatmap
Heatmap <- Heatmap(correlation_matrix_filtered_mat,
                   top_annotation = HeatmapAnnotation(
                     Stage = annotation_col$Stage,
                     col = list(Stage = annotation_colors$Stage),
                     show_legend = FALSE,
                     show_annotation_name = FALSE),
                   bottom_annotation = HeatmapAnnotation(
                     CellType = annotation_col$CellType,
                     col = list(CellType = annotation_colors$CellType),
                     show_annotation_name = FALSE),
                   col = colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
                   heatmap_legend_param = list(title = "Correlation",
                                               at = c(-0.5, 0, 0.5),
                                               labels = c("-0.5", "0", "0.5")),
                   cluster_rows = HCluster,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   row_split = 12,
                   column_split = factor(col_split, levels = unique(col_split))
)

# Show the heatmap
Heatmap <- draw(Heatmap)

# To export the heatmap as a .svg, use svg() function to open a connection to an SVG file. then draw the heatmap and close the connection using dev.off()
svg("HH8_HH9_SOX8_Coexpression_Heatmap.svg", width = 8, height = 6)
draw(Heatmap)
dev.off()

# Custom function to extract cluster gene IDs from a heatmap
Sox8_coexpression_clusters_HH8_HH9_0.1_12C <- Extract_HM_Clusters(Heatmap, correlation_matrix_filtered_mat)

# Adding a cluster label column to the original correlation_matrix
Sox8_coexpression_clusters_HH8_HH9_0.1_12C <- as.data.frame(Sox8_coexpression_clusters_HH8_HH9_0.1_12C[match(rownames(correlation_matrix_filtered_mat),Sox8_coexpression_clusters_HH8_HH9_0.1_12C[,1]), ])
correlation_matrix_filtered_df <- as.data.frame(correlation_matrix_filtered_mat)
correlation_matrix_filtered_df$Heatmap_Cluster <- c(Sox8_coexpression_clusters_HH8_HH9_0.1_12C[,2])

# Manually annotating Clusters based on co-expression at HH8 and HH9
 Cluster_annotations <- t(data.frame(
  Cluster_1 = c("Placodal", "Placodal"),
  Cluster_2 = c("Placodal", "Placodal"),
  Cluster_3 = c("NC", "Both"),
  Cluster_4 = c("NC", "NC"),
  Cluster_5 = c("NC", "NC"),
  Cluster_6 = c("Placodal", "Both"),
  Cluster_7 = c("NC", "Both"),
  Cluster_8 = c("NC", "NC"),
  Cluster_9 = c("Neither", "Placodal"),
  Cluster_10 = c("Placodal", "Placodal"),
  Cluster_11 = c("Both", "Both"),
  Cluster_12 = c("Placodal", "Neither")
 ))

# Using above dataframe to add cluster annotations to the correlation_matrix_dataframe
correlation_matrix_filtered_df <- cbind(correlation_matrix_filtered_df, Cluster_annotations[as.character(correlation_matrix_filtered_df[,"Heatmap_Cluster"]),])
colnames(correlation_matrix_filtered_df)[6:7] <- c("HH8","HH9")
correlation_matrix_filtered_df$Heatmap_Cluster <- as.numeric(sub(".*_(\\d+)$", "\\1", correlation_matrix_filtered_df$Heatmap_Cluster)) # Changing the cluster column to contain only numbers. This allows us to re-order the dataframe based on the order of the clusters
correlation_matrix_filtered_df <- correlation_matrix_filtered_df[order(correlation_matrix_filtered_df$Heatmap_Cluster),]

# Saving co-expression values together along with the cluster number and cell type annotations across stages
write.csv(correlation_matrix_filtered_df, "Heatmap_clusters/Sox8_coexpression_clusters_HH8_HH9_JASPAR_0.1_12C.csv")


#### PREPARE MEME FORMAT MOTIF LISTS BASED OFF CO-EXPRESSION ANALYSIS ABOVE ####

# Read in SOX8 co-expression heatmap clusters and extract TFs based on their co-expression patterns
HH8_HH9_Sox8_Coexpression_TFs <- read.csv("/data/Sox8_binding_partner_analysis/scRNAseq_objects/Heatmap_clusters/Sox8_coexpression_clusters_HH8_HH9_JASPAR_0.1_12C.csv")
colnames(HH8_HH9_Sox8_Coexpression_TFs)[1] <- "Gene"

# Instead of going off patterns from looking at the heatmap, I have selected TFs based on a Sox8 correlation coefficient of >0.1
HH8_Placodal_TFs <- HH8_HH9_Sox8_Coexpression_TFs[HH8_HH9_Sox8_Coexpression_TFs$HH8_Placodal > 0.1 & HH8_HH9_Sox8_Coexpression_TFs$HH8_NC < 0.1,]$Gene # >0.1 in placodes but not NC
HH8_NC_TFs <- HH8_HH9_Sox8_Coexpression_TFs[HH8_HH9_Sox8_Coexpression_TFs$HH8_NC > 0.1 & HH8_HH9_Sox8_Coexpression_TFs$HH8_Placodal < 0.1,]$Gene # >0.1 in NC but not placodes
HH8_NC_and_Placodal_TFs <- HH8_HH9_Sox8_Coexpression_TFs[HH8_HH9_Sox8_Coexpression_TFs$HH8_NC > 0.1 & HH8_HH9_Sox8_Coexpression_TFs$HH8_Placodal > 0.1,]$Gene # >0.1 in NC and Placodes
HH8_NC_and_Placodal_TFs <- append(HH8_NC_and_Placodal_TFs, "SOX8") # Adding Sox8 as a potential binding partner to this dataset
HH9_Placodal_TFs <- HH8_HH9_Sox8_Coexpression_TFs[HH8_HH9_Sox8_Coexpression_TFs$HH9_Placodal > 0.1 & HH8_HH9_Sox8_Coexpression_TFs$HH9_NC < 0.1,]$Gene
HH9_NC_TFs <- HH8_HH9_Sox8_Coexpression_TFs[HH8_HH9_Sox8_Coexpression_TFs$HH9_NC > 0.1 & HH8_HH9_Sox8_Coexpression_TFs$HH9_Placodal < 0.1,]$Gene
HH9_NC_and_Placodal_TFs <- HH8_HH9_Sox8_Coexpression_TFs[HH8_HH9_Sox8_Coexpression_TFs$HH9_NC > 0.1 & HH8_HH9_Sox8_Coexpression_TFs$HH9_Placodal > 0.1,]$Gene
HH9_NC_and_Placodal_TFs <- append(HH9_NC_and_Placodal_TFs, "SOX8") # Adding Sox8 as a potential binding partner to this dataset

# subset motif list based on motif lists of interest
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

write_meme(Sox8_motif, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/Sox8.txt", version = 5, ACGT_freqs)

write_meme(HH8_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_Placodal_motifs_0.1_min_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH8_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_NC_motifs_0.1_min_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH8_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH8_NC_and_Placodal_motifs_0.1_min_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH9_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH9_Placodal_motifs_0.1_min_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH9_NC_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH9_NC_motifs_0.1_min_PWM.txt", version = 5, ACGT_freqs)
write_meme(HH9_NC_and_Placodal_MotifList_vert, "/data/Sox8_binding_partner_analysis/meme_suite/MEME_motifs/HH9_NC_and_Placodal_motifs_0.1_min_PWM.txt", version = 5, ACGT_freqs)


