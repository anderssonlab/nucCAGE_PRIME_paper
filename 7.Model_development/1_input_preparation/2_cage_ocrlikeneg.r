#!/usr/bin/env Rscript
# get bw with provided genome info (to keep only std chr)

writeLines("\n #### 2 PREPARING TRAIN/TEST DATASET ####")

writeLines("\nImporting R libraries..")
suppressPackageStartupMessages({
  library(tools)
  library(argparse)
  library(tidyverse)
  library(dplyr)
  library(CAGEfightR)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(GenomicRanges)
})

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-i", "--input_dir",
                    help = "CAGE BigWig directory")
parser$add_argument("-r", "--ocr_rdata",
                    help = "DHS OCR RData file")
parser$add_argument("-m", "--design_matrix_file", default = "design.matrix.tsv",
                    help = "Design matrix of all CAGE data in .tsv format")
parser$add_argument("-n", "--name", default = "name", help = "Name")
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-d", "--distance", default = 200, type = "integer",
                    help = "Distance from the ATAC peak (bps)")
parser$add_argument("-t", "--test_chr", default = "chr2,chr3,chr4",
                    help = "Chromosome(s) that you would like to keep for the testing set (no space between , )") # nolint: line_length_linter.
parser$add_argument("-g", "--ucscgenome", default = "hg38",
                    help = "UCSC genome")
parser$add_argument("-s", "--species", default = "Homo sapiens",
                    help = "Species")

args <- parser$parse_args()

# General
dist <- args$distance
test_chr <- unlist(strsplit(args$test_chr, ","))

# Input and output
input_dir <- args$input_dir
ocr_rdata_file <- args$ocr_rdata
design_matrix_file <- args$design_matrix_file
output_dir <- args$output_dir
name <- args$name

ucsc_genome <- args$ucscgenome
species <- args$species



# 0. Helper functions
source("functions.r")

#' Extract CAGE Transcription Start Sites (CTSSs) from BigWig Files using CAGEfightR
get_ctss_from_bw <- function(dir_cage_bw,
                             design_matrix,
                             genome_info = NULL,
                             drop_chr = NULL,
                             keep_chr = NULL) {

  # Create BigWigFileList objects for plus and minus strands
  bwplus <- rtracklayer::BigWigFileList(file.path(dir_cage_bw,
                                                  design_matrix$BigWigPlus))
  bwminus <- rtracklayer::BigWigFileList(file.path(dir_cage_bw,
                                                   design_matrix$BigWigMinus))
  names(bwplus) <- names(bwminus) <- rownames(design_matrix)

  # Handle genome_info, drop_chr, and keep_chr
  gn <- genome_info
  if (!is.null(genome_info)) {
    if (!is.null(drop_chr)) {
      gn <- GenomeInfoDb::dropSeqlevels(gn, drop_chr)
    }
    if (!is.null(keep_chr)) {
      gn <- GenomeInfoDb::keepSeqlevels(gn, keep_chr)
    }
  } else if (!is.null(drop_chr) || !is.null(keep_chr)) {
    stop("genome_info must be provided when using drop_chr or keep_chr.")
  }

  # Quantify CTSSs, calculate pooled counts and TPM
  orig_ctss <- CAGEfightR::quantifyCTSSs(plusStrand = bwplus,
                                         minusStrand = bwminus,
                                         design = design_matrix,
                                         genome = gn)
  orig_ctss <- CAGEfightR::calcPooled(orig_ctss, inputAssay = "counts")
  orig_ctss <- CAGEfightR::calcTPM(orig_ctss,
                                   inputAssay = "counts",
                                   outputAssay = "TPM",
                                   totalTags = NULL,
                                   outputColumn = "totalTags")
  return(orig_ctss)
}

# Character of negative
# 1. Singleton was removed before calling TC
# 2. TC won't overlap with the positive ocr
# 3. Extend DIST from THICK of TC
# It means that TC won't overlap with ocr, but extended TC might
# = we allow the neg to overlap with pos
# Note that TCs come from all selected CAGE.

#' Generate OCR-like Negative Samples
ocrlike_for_negative <- function(orig_ctss, ocr_gr, dis) {
  supp_ctss <- CAGEfightR::subsetBySupport(orig_ctss,
                                           inputAssay = "counts",
                                           outputColumn = "support",
                                           unexpressed = 1,
                                           minSamples = 0)
  supp_ctss <- CAGEfightR::calcPooled(supp_ctss,
                                      inputAssay = "counts")
  tcs <- CAGEfightR::clusterUnidirectionally(supp_ctss)
  tcs_neg <- IRanges::subsetByOverlaps(tcs,
                                       ocr_gr,
                                       invert = TRUE)

  # Extend the negative (ocr-like) positions from the center of the thick position # nolint: line_length_linter.
  tcs_neg <- extend_from_center_thick_gr(tcs_neg, dis, keep_same_length = TRUE)
  mcols(tcs_neg) <- mcols(tcs_neg)[, "thick", drop = FALSE]

  # update row names to remove strand information
  strand(tcs_neg) <- S4Vectors::Rle(rep("*", length(tcs_neg)))

  # Remove duplicated genomic ranges
  duplicated_indices <- GenomicRanges::duplicated(tcs_neg)

  # Subset the GRanges object to keep only unique ranges
  print("Removing duplicated genomic ranges")
  unique_tcs_neg <- tcs_neg[!duplicated_indices]
  names(unique_tcs_neg) <- paste0(seqnames(unique_tcs_neg), ":",
                                  start(unique_tcs_neg), "-",
                                  end(unique_tcs_neg), ";",
                                  strand(unique_tcs_neg))

  # Report
  len_neg <- length(unique_tcs_neg)
  print(paste0("Number of ocr-like positions: ", len_neg))
  print(paste0("Ocr-like positions that (any) overlaps with the ocr: ", length(IRanges::subsetByOverlaps(unique_tcs_neg, ocr_gr, type = c("any")))))  # nolint: line_length_linter.
  print(paste0("Sanity check to make sure that no exact overlap of negative (ocr-like) to positive (ocr) (expect to be 0): ", length(IRanges::subsetByOverlaps(unique_tcs_neg, ocr_gr, type = c("equal"))))) # nolint: line_length_linter.

  return(unique_tcs_neg)

}



# 1. Load RData
load(ocr_rdata_file)



# 2. Read in CAGE training and testing data, based on CAGEFightR,
# no filtering/isubsetBySupport
writeLines("\nReading in data")
design_matrix <- read.table(design_matrix_file,
                            header = TRUE,
                            sep = "\t",
                            row.names = "Name")
genome_info <- GenomeInfoDb::keepStandardChromosomes(rtracklayer::SeqinfoForUCSCGenome(ucsc_genome), # nolint: line_length_linter.
                                                     species = species)
genome_info <- GenomeInfoDb::dropSeqlevels(genome_info, "chrM")

tr_orig_ctss <- get_ctss_from_bw(input_dir,
                                 design_matrix,
                                 genome_info,
                                 drop_chr = test_chr)
te_orig_ctss <- get_ctss_from_bw(input_dir,
                                 design_matrix,
                                 genome_info,
                                 keep_chr = test_chr)



# 3. Create negative from TC (ocr-like)

# Character of negative
# 1. Singleton was removed before calling TC
# 2. TC won't overlap with the positive ocr
# 3. Extend DIST from THICK of TC
# It means that TC won't overlap with ocr, but extended TC might
# = we allow the neg to overlap with pos
# Note that TCs come from all selected CAGE.

writeLines("\nCreating the ocr-like from TCs of training data..")
ocrlike_neg_train_gr <- ocrlike_for_negative(tr_orig_ctss, ocr_train_gr, dist)
writeLines("\nCreating the ocr-like from TCs of testing data..")
ocrlike_neg_test_gr <- ocrlike_for_negative(te_orig_ctss, ocr_test_gr, dist)



# 4. Save
writeLines("\nSaving..")
save(list = c("ocr_train_gr",
              "ocr_test_gr",
              "ocrlike_neg_train_gr",
              "ocrlike_neg_test_gr"),
     file = paste0("2_", name, "_ocr_ocrlikeneg.RData"))
save(list = c("tr_orig_ctss", "te_orig_ctss"),
     file = paste0("2_", name, "_orig_ctss.RData"))
