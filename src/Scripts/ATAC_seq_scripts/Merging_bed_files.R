library(GenomicRanges)
library(rtracklayer)

# Define file paths
bed_files <- c("/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/HH9_NC_Sox8_Sox4_2bp.bed", 
               "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/HH9_NC_Sox8_Sox4_3bp.bed", 
               "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/HH9_NC_Sox8_Sox4_4bp.bed")  # Replace with actual file paths

# Read BED files as GRanges and standardize seqlevels
gr_list <- lapply(bed_files, function(f) {
  gr <- import(f, format = "BED")  # Import as GRanges
  seqlevelsStyle(gr) <- "UCSC"  # Standardize chromosome names (e.g., "chr1", "chr2")
  return(gr)
})

# Ensure all GRanges have the same seqlevels
common_seqlevels <- Reduce(intersect, lapply(gr_list, seqlevels))
gr_list <- lapply(gr_list, function(gr) keepSeqlevels(gr, common_seqlevels, pruning.mode="coarse"))

# Merge overlapping regions and remove duplicates
merged_gr <- reduce(do.call(c, gr_list))

# Write the final merged BED file
export(merged_gr, "merged_output.bed", format = "BED")