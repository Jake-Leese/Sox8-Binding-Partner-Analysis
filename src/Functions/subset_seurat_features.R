# Function to subset Seurat objects to include either JASPAR2024 or AnimalDB TFs

subset_seurat_features <- function(Seurat_obj, DB) {
  
  if (DB == "JASPAR2024") {
    
    Jaspar2024 <- JASPAR2024()
    sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(Jaspar2024))

    motifList <- getMatrixSet(x = sq24, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

    # rename each motif from their ID to their TF name
    name_vector <- c()
    for (i in 1:length(motifList)){
      name <- name(motifList[[i]])
      name_vector <- c(name_vector, name)
    }
    names(motifList) <- name_vector
    
    # Subsetting Seurat object to include JASPAR TFs only
    Seurat_obj_TFs <- DietSeurat(Seurat_obj, features = name_vector)
    
    return(Seurat_obj_TFs)
    
  } else if (DB == "AnimalTFDB") {
    
    GalGal_TFs_df <- read.table("/data/Sox8_binding_partner_analysis/Database_objects/Gallus_gallus_TF.txt", sep = "\t", header = TRUE)
    GalGal_TFs <- GalGal_TFs_df$Symbol
    
    Seurat_obj_TFs <- DietSeurat(Seurat_obj, features = GalGal_TFs)
    
    return(Seurat_obj_TFs)
  }
}
