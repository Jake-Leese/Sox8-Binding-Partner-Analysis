# Script to:
# 1. read in .bed files from e.g. SpaMo outputs
# 2. Match .bed ranges to ArchR peaks
# 3. Check overlaps with peak2gene links to identify downstream targets
# 3. Filter downstream targets by those with placodal expression (based on scRNA-seq data)

# Written using alexthiery-schelper-archr_dev_macs2-schelper-0.3.5.img 

.libPaths("/R/libs/AT_ArcR_macs2")
setwd('/data/Sox8_binding_partner_analysis/scATACseq_objects')

library(ArchR)
library(tidyverse)
library(BSgenome.Ggallus.UCSC.galGal6)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(universalmotif)
library(PANTHER.db)
library(biomaRt)


source("/data/Sox8_binding_partner_analysis/src/Functions/convert_coordinates_to_GRanges.R")

bed_file_path <- "/data/Sox8_binding_partner_analysis/enhancer_annotation/bed_files/PPR/Reproducible_open/"

# Load FullData ArchR project with peak2gene links
FullData <- loadArchRProject(path = '/data/scATAC_paper/FullData_Save-ArchR')

# Function to retrieve downstream genes to ArchR peaks based on p2g links. The outputted object is a dataframe with gene name, ensembl_gene_id, hgnc_symbol, go_id, go_term_name
match_peaks_to_ds_genes <- function(BED, ArchRProj, p2g_corr_cutoff, p2g_FDR_cutoff) {

##############################################################################################################################
############################################### READ IN PEAKS ################################################################

# Function to add "chr" prefix to seqnames (this is necessary to check overlaps with peaksets from ArchR)
add_chr_prefix <- function(gr) {
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  
  seqlevels(gr) <- ifelse(grepl("^chr", seqlevels(gr)), seqlevels(gr), paste0("chr", seqlevels(gr)))
  seqnames(gr) <- factor(seqnames(gr), levels = seqlevels(gr))
  
  return(gr)
}

# Import .bed file as GRanges and append the seqnames with "chr"
BED.gr <- add_chr_prefix(import(BED, format = "BED"))


##############################################################################################################################
########################################## READ IN PEAK2GENE LINKS ###########################################################

# Get peak2GeneLinks use corCutOff = 0 to retrieve any positively correlated genes (This is what was used for GRNi. Less stringent) 0.45 is ArchR default
p2g <- getPeak2GeneLinks(ArchRProj, corCutOff = p2g_corr_cutoff, FDRCutOff = p2g_FDR_cutoff, returnLoops = FALSE)

# Convert peak2gene links dataframe to a GRanges object for further cross comparisons
p2g_gr <- metadata(p2g)$peakSet[p2g$idxATAC]
mcols(p2g_gr)$linked_gene <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]


##############################################################################################################################
############################################### CHECK OVERLAP ################################################################

# Matching original ArchR peakset to imported bed peaks (based on overlap)
peakset <- getPeakSet(ArchRProj)
peakset_overlaps <- findOverlaps(peakset, BED.gr)
archr_bed_peakset <- peakset[unique(queryHits(peakset_overlaps))]

# Get correlated linked genes 
linked_genes <- unique(p2g_gr[p2g_gr %in% archr_bed_peakset]$linked_gene)


##############################################################################################################################
########################### Converting gene names into supported format for pantherDB ########################################

ggallus_ensembl <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl", version = 106, mirror = "useast") # extracting ensembl for galgal6

# Query the database for UniProt IDs from downstream genes
results <- getBM(
  attributes = c("external_gene_name","ensembl_gene_id", "hgnc_symbol", "uniprot_gn_id", "go_id", "name_1006"),
  filters = "external_gene_name",
  values = linked_genes,
  mart = ggallus_ensembl
)

return(results)
}


HH9_IKZF2_ds_genes_up <- match_peaks_to_ds_genes(
  BED = paste0(bed_file_path, "HH9_PPR_rep_SOX8_IKZF2_0bp_up.bed"),
  ArchRProj = FullData,
  p2g_corr_cutoff = 0.2,
  p2g_FDR_cutoff = 1e-04)


write_lines(unique(HH9_IKZF2_ds_genes_up$external_gene_name), "/data/Sox8_binding_partner_analysis/enhancer_annotation/p2g_outputs/HH9_IKZF2_ds_genes_up_0.2_corr.txt")



