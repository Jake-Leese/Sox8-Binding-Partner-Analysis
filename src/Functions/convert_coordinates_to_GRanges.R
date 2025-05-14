# Function that takes a character vector of genome co-ordinates in "chr-start-end" format, and returns each co-ordinate as a range in a GRanges object

convert_coordinates_to_GRanges <- function(coordinates){

DF <- data.frame(coordinates)                    # Create a dataframe, DF, that separates seqnames (chrx) from genomic ranges.
DF <- DF %>%                                    # Issue with this, is that it separates out values by the presence of a "-"
  separate(                                     # Ranges column therefore only contains the first value of the IRange before the "-"
    coordinates,
    into = c("seqnames", "ranges"),
    sep = "-",
    extra = "drop")
seqnames <- as.vector(DF$seqnames) # Save seqnames for GRanges

# EXTRACT CHARACTERS AFTER FIRST DASH
# Initialize an empty vector to store the results
  result_vector <- character(length(coordinates))
  
  # Loop through each element in the input vector
  for (i in seq_along(coordinates)) {
    # Find the position of the first "-"
    first_dash_position <- regexpr("-", coordinates[i])
    
    # Check if a dash is found
    if (first_dash_position != -1) {
      # Extract the substring starting from the position after the first "-"
      result_vector[i] <- substr(coordinates[i], first_dash_position + 1, nchar(coordinates[i]))
    } else {
      # If no dash is found, retain the entire string
      result_vector[i] <- coordinates[i]
    }
  }


# Create a new GRanges object using the seqnames and ranges
GR <- GRanges(seqnames = seqnames, ranges = result_vector)

return(GR)
}
