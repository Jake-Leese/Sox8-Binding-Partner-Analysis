# Function to extract count data from Seurat object. 
# Comes with options for stage, sample, assay (RNA vs sct), and data type (counts vs logNormalised)

Extract_count_data <- function(Seurat_obj, Stage, Cell_type, Sample, Assay, data_type) {
  
  # SUBSET CELLS
  Seurat_processed <- Seurat_obj[,Seurat_obj$ectoderm_type %in% Cell_type]
  Seurat_processed <- Seurat_processed[,Seurat_processed$stage %in% Stage]
  
  if (!missing(Sample)) {
    # Action taken when sample is specified
    Seurat_processed <- Seurat_processed[,Seurat_processed$orig.ident %in% Sample]
  } else {
    NULL # do nothing
  }
  
  # EXTRACT COUNTS
  Seurat_processed <- JoinLayers(Seurat_processed, assay = "RNA") # Joining RNA count layers
  
  if (missing(data_type)) {
    Seurat_processed <- Seurat_processed[[Assay]]$counts
    return(Seurat_processed)
  } else {
    Seurat_processed <- GetAssayData(Seurat_processed, assay = Assay, slot = data_type)
    return(Seurat_processed)
  }
}

