# Function to convert a JASPAR2024 PWMatrix into a MEME minimal motif format
# Function to write a MEME simple format motif file from JASPAR
Jaspar2MEME <- function(obj, file, A, C, G, Th) {
  # Extract the slot names of the S4 object
  slot_names <- slotNames(obj)
  
  # Open a connection to the text file
  con <- file(file, "w")
  
  writeLines("MEME version 4", con)
  writeLines("", con)
  writeLines("ALPHABET= ACTG", con)
  writeLines("", con)
  writeLines("strands: + -", con)
  writeLines("", con)
  writeLines("Background letter frequencies", con)
  char_vector <- c("A", "C", "G", "T")
  value_vector <- c(A, C, G, Th)
  output_string <- ""
  for (i in seq_along(char_vector)) {
    if (i > 1) {
      output_string <- paste(output_string, char_vector[i], value_vector[i], sep = " ")
    } else {
      output_string <- paste(output_string, char_vector[i], value_vector[i], sep = "")
    }
  }
  writeLines(output_string, con)
  writeLines("", con)
  Motif_ID <- c("MOTIF", paste0(obj@ID), paste0(obj@name))
  Motif_ID <- paste(Motif_ID, collapse = " ")
  writeLines(Motif_ID, con)
  for (i in seq_len(ncol(obj@profileMatrix))) {
    writeLines(paste(paste0(obj@profileMatrix[,i]), collapse = " "), con)
  }
  
  # Close the file connection
  close(con)
}