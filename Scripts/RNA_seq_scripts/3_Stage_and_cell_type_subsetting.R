# Subsetting gene list in original HHall_ectoderm object to include just TFs. 
# Further subsetting cells by stage and ectoderm types. And saving as seurat rds objects. 


setwd("/data/scRNAseq_objects")
.libPaths("/R/libs/R_scvi_integration_v3.4")

library(Seurat)
library(usethis)
library(devtools)

# Load and check the Seurat object
HHall_ectoderm <- readRDS("HHall_ectoderm_2024-06-24")

# Importing GalGal_TF list from AnimalTFDB4
GalGal_TFs_df <- read.table("Gallus_gallus_TF.txt", sep = "\t", header = TRUE)
GalGal_TFs <- GalGal_TFs_df$Symbol

# Subsetting genes (features) to include GalGal TFs only
HHall_ectoderm_TFs <- DietSeurat(HHall_ectoderm, features = GalGal_TFs)
SaveSeuratRds(HHall_ectoderm_TFs, file = "HHall_ectoderm_TFs")


# Extract all NC cells, subset by stage (HH7 - HH16), and save objects to disk
HHall_NC_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$ectoderm_type %in% "NC"]

HH7_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH7"]
SaveSeuratRds(HH7_NC_TFs, file = "HH7_NC_TFs")

HH8_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH8"]
SaveSeuratRds(HH8_NC_TFs, file = "HH8_NC_TFs")

HH9_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH9"]
SaveSeuratRds(HH9_NC_TFs, file = "HH9_NC_TFs")

HH12_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH12"]
SaveSeuratRds(HH12_NC_TFs, file = "HH12_NC_TFs")

HH14_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH14"]
SaveSeuratRds(HH14_NC_TFs, file = "HH14_NC_TFs")

HH16_NC_TFs <- HHall_NC_TFs[,HHall_NC_TFs$stage %in% "HH16"]
SaveSeuratRds(HH16_NC_TFs, file = "HH16_NC_TFs")


# Extract all placodal cells, subset by stage (HH7 - HH16), and save objects to disk
HHall_Placodal_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$ectoderm_type %in% "placode"]

HH7_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH7"]
SaveSeuratRds(HH7_Placodal_TFs, file = "HH7_Placodal_TFs")

HH8_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH8"]
SaveSeuratRds(HH8_Placodal_TFs, file = "HH8_Placodal_TFs")

HH9_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH9"]
SaveSeuratRds(HH9_Placodal_TFs, file = "HH9_Placodal_TFs")

HH12_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH12"]
SaveSeuratRds(HH12_Placodal_TFs, file = "HH12_Placodal_TFs")

HH14_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH14"]
SaveSeuratRds(HH14_Placodal_TFs, file = "HH14_Placodal_TFs")

HH16_Placodal_TFs <- HHall_Placodal_TFs[,HHall_Placodal_TFs$stage %in% "HH16"]
SaveSeuratRds(HH16_Placodal_TFs, file = "HH16_Placodal_TFs")


# Extract NC & Placodal cells, subset by stage (HH7 - HH16), and save objects to disk
HHall_NC_and_Placodal_TFs <- HHall_ectoderm_TFs[,HHall_ectoderm_TFs$ectoderm_type %in% c("NC", "placode")]

HH7_NC_and_Placodal_TFs <- HHall_NC_and_Placodal_TFs[,HHall_NC_and_Placodal_TFs$stage %in% "HH7"]
SaveSeuratRds(HH7_NC_and_Placodal_TFs, file = "HH7_NC_and_Placodal_TFs")

HH8_NC_and_Placodal_TFs <- HHall_NC_and_Placodal_TFs[,HHall_NC_and_Placodal_TFs$stage %in% "HH8"]
SaveSeuratRds(HH8_NC_and_Placodal_TFs, file = "HH8_NC_and_Placodal_TFs")

HH9_NC_and_Placodal_TFs <- HHall_NC_and_Placodal_TFs[,HHall_NC_and_Placodal_TFs$stage %in% "HH9"]
SaveSeuratRds(HH9_NC_and_Placodal_TFs, file = "HH9_NC_and_Placodal_TFs")

HH12_NC_and_Placodal_TFs <- HHall_NC_and_Placodal_TFs[,HHall_NC_and_Placodal_TFs$stage %in% "HH12"]
SaveSeuratRds(HH12_NC_and_Placodal_TFs, file = "HH12_NC_and_Placodal_TFs")

HH14_NC_and_Placodal_TFs <- HHall_NC_and_Placodal_TFs[,HHall_NC_and_Placodal_TFs$stage %in% "HH14"]
SaveSeuratRds(HH14_NC_and_Placodal_TFs, file = "HH14_NC_and_Placodal_TFs")

HH16_NC_and_Placodal_TFs <- HHall_NC_and_Placodal_TFs[,HHall_NC_and_Placodal_TFs$stage %in% "HH16"]
SaveSeuratRds(HH16_NC_and_Placodal_TFs, file = "HH16_NC_and_Placodal_TFs")

# Checked each subset and sample numbers add up in every case

