#!/usr/bin/env Rscript

writeLines("\n #### 3 DATA AUGMENTATION ####")

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

set.seed(214)

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = "./",
                    help = "OCR and OCR-like RData file")
parser$add_argument("-n", "--name", default = "name", help = "Name")
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "output directory")
parser$add_argument("-d", "--distance", default = 200, type = "integer",
                    help = "distance from the DHS peak (bps)")

args <- parser$parse_args()

# general
dist <- args$distance

# input and output
ocr_rdata_file <- args$ocr_rdata
output_dir <- args$output_dir
name <- args$name



# 0 helper functions
source("functions.r")

#' Generate Random Integers with Normal Distribution
rnorm_int <- function(number) {
  ch <- 0
  ls <- c()

  while (ch < number) {
    cntrl_num <- number - ch
    nn <- round(stats::rnorm(cntrl_num, mean = 0, sd = 10), 0)
    nn <- nn[nn >= -30 & nn <= 30]
    ls <- c(ls, nn)
    ch <- length(ls)
  }

  return(ls[1:number])
}

#' Generate Random Numbers with Uniform and Normal Distributions
get_random_number <- function(n) {
  unif_int <- round(runif(n = n, min = -30, max = 30), 0)

  norm_int <- rnorm_int(n)

  return(list("unif" = unif_int, "norm" = norm_int))
}

#' Data Augmentation with Uniform Distribution
augmentation_unif <- function(ocr) {
  # Assert that the 'thick' column is present and is an IRanges object
  assertthat::assert_that("thick" %in% names(S4Vectors::mcols(ocr)),
                          msg = "The 'thick' column is not present in the GRanges object.") # nolint: line_length_linter.
  assertthat::assert_that(is(ocr$thick, "IRanges"),
                          msg = "The 'thick' column is not an IRanges object.")

  num <- length(ocr)
  set <- get_random_number(num)

  # ocr_unif
  ocr_unif <- ocr
  start(ocr_unif) <- start(ocr_unif) + set$unif
  end(ocr_unif) <- end(ocr_unif) + set$unif
  ocr_unif$thick <- IRanges::shift(ocr_unif$thick, shift = set$unif)
  ocr_unif <- IRanges::subsetByOverlaps(ocr_unif,
                                        ocr,
                                        type = "equal",
                                        invert = TRUE)

  return(ocr_unif)
}

#' Data Augmentation with Normal Distribution
augmentation_norm <- function(ocr) {
  # Assert that the 'thick' column is present and is an IRanges object
  assertthat::assert_that("thick" %in% names(S4Vectors::mcols(ocr)),
                          msg = "The 'thick' column is not present in the GRanges object.") # nolint: line_length_linter.
  assertthat::assert_that(is(ocr$thick, "IRanges"),
                          msg = "The 'thick' column is not an IRanges object.")

  num <- length(ocr)
  set <- get_random_number(num)

  # ocr_norm
  ocr_norm <- ocr
  start(ocr_norm) <- start(ocr_norm) + set$norm
  end(ocr_norm) <- end(ocr_norm) + set$norm
  ocr_norm$thick <- IRanges::shift(ocr_norm$thick, shift = set$norm)
  ocr_norm <- IRanges::subsetByOverlaps(ocr_norm,
                                        ocr,
                                        type = "equal",
                                        invert = TRUE)

  return(ocr_norm)
}

nearpos_negative_front <- function(ocr_gr, dist, core) {
  # Shift core data and 'thick' column by -(200 + 75 + 1)
  nearpos_front_gr <- IRanges::shift(ocr_gr, shift = -(dist + core + 1))
  nearpos_front_gr$thick <- IRanges::shift(nearpos_front_gr$thick,
                                           shift = -(dist + core + 1))

  return(nearpos_front_gr)
}

nearpos_negative_back <- function(ocr_gr, dist, core) {
  # Shift core data and 'thick' column by (200 + 75 + 1)
  nearpos_back_gr <- IRanges::shift(ocr_gr, shift = (dist + core + 1))
  nearpos_back_gr$thick <- IRanges::shift(nearpos_back_gr$thick,
                                          shift = (dist + core + 1))

  return(nearpos_back_gr)
}

get_core_gr <- function(gr, dist, core) {
  core_gr <- GRanges(seqnames = seqnames(gr),
                     ranges = IRanges(start = start(gr) + (dist - core),
                                      end = end(gr) - (dist - core)))
  core_gr$thick <- gr$thick
  names(core_gr) <- names(gr)
  return(core_gr)
}

finalize_nearpos_neg <- function(npn_front_gr,
                                 npn_back_gr,
                                 ocr_core_gr,
                                 ocr_aug_core_gr,
                                 dist,
                                 core) {

  # npn = near-positive negatives

  # Remove npn won't overlaps with positive core
  npn_front_gr <- IRanges::subsetByOverlaps(npn_front_gr,
                                            ocr_core_gr,
                                            invert = TRUE)
  npn_back_gr <- IRanges::subsetByOverlaps(npn_back_gr,
                                           ocr_core_gr,
                                           invert = TRUE)

  # "CORE" np-pos neg won't overlap with core pos aug

  #### "CORE" npn
  npn_front_core <- get_core_gr(npn_front_gr, dist, core)
  npn_back_core <- get_core_gr(npn_back_gr, dist, core)

  #### won't overlap with positive augmented core
  npn_front_core <- IRanges::subsetByOverlaps(npn_front_core,
                                              ocr_aug_core_gr,
                                              invert = TRUE)
  npn_back_core <- IRanges::subsetByOverlaps(npn_back_core,
                                             ocr_aug_core_gr,
                                             invert = TRUE)

  #### resize gr
  npn_front <- expand_core_gr(npn_front_core, dist, core)
  npn_back <- expand_core_gr(npn_back_core, dist, core)

  # Combine front and back cores
  combined_npn_gr <- c(npn_front, npn_back)
  combined_npn_gr <- sort(combined_npn_gr)

  names(combined_npn_gr) <- paste0(seqnames(combined_npn_gr), ":",
                                   start(combined_npn_gr), "-",
                                   end(combined_npn_gr), ";",
                                   strand(combined_npn_gr))

  return(combined_npn_gr)
}

expand_core_gr <- function(core_gr, dist, core) {
  gr <- GRanges(seqnames = seqnames(core_gr),
                ranges = IRanges(start = start(core_gr) - (dist - core),
                                 end = end(core_gr) + (dist - core)))
  gr$thick <- core_gr$thick
  names(gr) <- names(core_gr)

  return(gr)
}



# 1 load Rdata
load(ocr_rdata_file)



# 2 augment the positive OCR based on
# normal distribution prob within +- 30 bps from center
# they were check to make sure that
# the augment one not exact overlap with original
writeLines("\naugmenting the training- positive- ocr based on normal distribution prob within +- 30 bps from center..") # nolint: line_length_linter.
ocr_train_aug_gr <- augmentation_norm(ocr_train_gr)
print(paste0("train positive: ", length(ocr_train_gr)))
print(paste0("augmented train positive: ", length(ocr_train_aug_gr)))

writeLines("\naugmenting the testing- positive- ocr based on normal distribution prob within +- 30 bps from center..") # nolint: line_length_linter.
ocr_test_aug_gr <- augmentation_norm(ocr_test_gr)
print(paste0("test positive: ", length(ocr_test_gr)))
print(paste0("augmented test positive: ", length(ocr_test_aug_gr)))



# 3 augment the negative OCR based on
# random distribution prob within +- 30 bps from center
# they were check to make sure tha
# the augment one not exact overlap with original
writeLines("\naugmenting the trianing- negative- ocr based on random distribution prob within +- 30 bps from center..") # nolint: line_length_linter.
ocrlike_neg_train_aug_gr <- augmentation_unif(ocrlike_neg_train_gr)
print(paste0("train negative: ", length(ocrlike_neg_train_gr)))
print(paste0("augmented train negative (before filtering): ",
             length(ocrlike_neg_train_aug_gr)))

writeLines("\naugmenting the testing- negative- ocr based on random distribution prob within +- 30 bps from center..") # nolint: line_length_linter.
ocrlike_neg_test_aug_gr <- augmentation_unif(ocrlike_neg_test_gr)
print(paste0("test negative: ", length(ocrlike_neg_test_gr)))
print(paste0("augmented test negative (before filtering): ",
             length(ocrlike_neg_test_aug_gr)))



# 4 subset the training negative to training positive (orig + aug)
# we allow the negative to be "50% or less" overlap with the positive
ocrlike_neg_train_aug_gr <- IRanges::subsetByOverlaps(ocrlike_neg_train_aug_gr,
                                                      c(ocr_train_aug_gr,
                                                        ocr_train_gr),
                                                      minoverlap = dist,
                                                      invert = TRUE)
print(paste0("\naugmented train negative (remove more than 50% overlap with pos): ", # nolint: line_length_linter.
             length(ocrlike_neg_train_aug_gr)))

ocrlike_neg_test_aug_gr <- IRanges::subsetByOverlaps(ocrlike_neg_test_aug_gr,
                                                     c(ocr_test_aug_gr,
                                                       ocr_test_gr),
                                                     minoverlap = dist,
                                                     invert = TRUE)
print(paste0("\naugmented train negative (remove more than 50% overlap with pos): ", # nolint: line_length_linter.
             length(ocrlike_neg_test_aug_gr)))




# 5 set up near-positve negative,
# core data +/- 75 bps

core <- 75

npn_front_tr <- nearpos_negative_front(ocr_train_gr, dist, core)
npn_back_tr <- nearpos_negative_back(ocr_train_gr, dist, core)

npn_front_te <- nearpos_negative_front(ocr_test_gr, dist, core)
npn_back_te <- nearpos_negative_back(ocr_test_gr, dist, core)

ocr_train_core <- get_core_gr(ocr_train_gr, dist, core)
ocr_train_aug_core <- get_core_gr(ocr_train_aug_gr, dist, core)

ocr_test_core <- get_core_gr(ocr_test_gr, dist, core)
ocr_test_aug_core <- get_core_gr(ocr_test_aug_gr, dist, core)

nearpos_neg_train_aug_gr <- finalize_nearpos_neg(npn_front_tr,
                                                 npn_back_tr,
                                                 ocr_train_core,
                                                 ocr_train_aug_core,
                                                 dist,
                                                 core)

nearpos_neg_test_aug_gr <- finalize_nearpos_neg(npn_front_te,
                                                npn_back_te,
                                                ocr_test_core,
                                                ocr_test_aug_core,
                                                dist,
                                                core)



# 6 save
writeLines("\nSaving..")
save(list = c("ocr_train_gr",
              "ocr_train_aug_gr",
              "ocr_test_gr",
              "ocr_test_aug_gr",
              "ocrlike_neg_train_gr",
              "ocrlike_neg_train_aug_gr",
              "ocrlike_neg_test_gr",
              "ocrlike_neg_test_aug_gr",
              "nearpos_neg_train_aug_gr",
              "nearpos_neg_test_aug_gr"),
     file = paste0("3_", name, "_augmentation.RData"))
