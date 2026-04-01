#!/usr/bin/env Rscript

writeLines("\n #### 5 VALIDATION AND SELECTED NEGATIVE ####")

writeLines("\nImporting R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(dplyr)
  library(CAGEfightR)
  library(AnnotationDbi)
  library(GenomicFeatures)

  library(rtracklayer)
  library(tibble)
  library(readr)
  library(forcats)
  library(tidyr)
  library(GenomicRanges)
  library(reshape2)
})

set.seed(214)

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = "./",
                    help = "OCR and OCR-like RData file ")
parser$add_argument("-c", "--cage_tss_rdata", default = ".",
                    help = "CAGE TSS rData file")
parser$add_argument("-n", "--name", default = "name",
                    help = "Name")
parser$add_argument("-v", "--validation_ratio", default = 0.2,
                    help = "validation_ratio")

args <- parser$parse_args()

# input and output
ocr_rdata_file <- args$ocr_rdata
tss_rdata_file <- args$cage_tss_rdata
name <- args$name
validation_ratio <- args$validation_ratio



#0 helper functions
source("functions.r")

select_ocr_on_condition_each_lib <- function(cage_rse,
                                             ocr_gr,
                                             coln,
                                             min_score_for_onebase_ocr = 4) {
  cage_gr <- cast_rse_to_granges(cage_rse, assay = "counts", coln_assay = coln)
  cage_filtered_gr <- cage_gr[cage_gr$score > 0]
  expr <- GenomicRanges::countOverlaps(ocr_gr, cage_filtered_gr)
  ocr_df <- data.frame(ocr_gr)
  ocr_df$num_overlaps <- expr

  sure_ocr_df <- dplyr::filter(ocr_df, num_overlaps > 1)

  onebase_ocr_df <- dplyr::filter(ocr_df, num_overlaps == 1)
  onebase_ocr_gr <- GenomicRanges::GRanges(onebase_ocr_df)
  onebase_ovl <- GenomicRanges::findOverlaps(onebase_ocr_gr, cage_filtered_gr)

  cage_hit_gr <- cage_filtered_gr[S4Vectors::subjectHits(onebase_ovl), ]
  onebase_ocr_df$onebase_score <- cage_hit_gr$score

  slt_onebase_ocr_df <- dplyr::filter(onebase_ocr_df, onebase_score >= min_score_for_onebase_ocr)
  slt_onebase_ocr_df <- dplyr::select(slt_onebase_ocr_df, -onebase_score)

  slt_ocr_df <- rbind(sure_ocr_df, slt_onebase_ocr_df)
  slt_ocr_df <- dplyr::arrange(slt_ocr_df, seqnames, start, end)

  return(GenomicRanges::GRanges(slt_ocr_df))
}

#' Select OCRs Based on Specific Conditions
select_ocr_on_condition <- function(cage_rse, ocr_gr) {
  slt_ocr_grl <- GenomicRanges::GRangesList()
  for (i in 1:ncol(cage_rse)) {
    slt_ocr_grl[[i]] <- select_ocr_on_condition_each_lib(cage_rse, ocr_gr, i)
  }
  return(slt_ocr_grl)
}

#' Subset Positive Annotations
subset_annotation_positive <- function(pos_sltocr_grl) {
  grl <- GenomicRanges::GRangesList()
  for (i in seq_along(pos_sltocr_grl)) {
    # Assert that txType is present in the GRanges object
    assertthat::assert_that("txType" %in% names(GenomicRanges::mcols(pos_sltocr_grl[[i]])),
                            msg = paste("txType is not present in the GRanges object at index", i))

    grl[[i]] <- S4Vectors::subset(pos_sltocr_grl[[i]], !(txType %in% c('exon', 'CDS', 'threeUTR')))
  }
  return(grl)
}

#' Select Negative Samples to Match Positive Samples
select_negative_to_match_positive <- function(neg_sltocr_grl, match_pos_grl) {
  grl <- GenomicRanges::GRangesList()

  for (i in seq_along(neg_sltocr_grl)) {
    new_df <- data.frame(neg_sltocr_grl[[i]])
    num_pos <- length(match_pos_grl[[i]])
    neg_ratio <- num_pos / nrow(new_df)

    split_data <- new_df %>%
      dplyr::group_by(seqnames, txType) %>%
      dplyr::group_split() %>%
      purrr::map(function(group_data) {
        sample_split <- dplyr::sample_frac(group_data, neg_ratio)
        list(negative = sample_split)
      })

    split_neg <- purrr::map(split_data, function(split_element) split_element$negative) %>% 
      dplyr::bind_rows()

    grl[[i]] <- GenomicRanges::GRanges(split_neg)
  }
  return(grl)
}

#' Split Data into Training and Validation Sets
split_train_valid <- function(df, ratio_valid = 0.2) {
  # split based on chr and txType
  split_data <- df %>%
    dplyr::group_by(seqnames, txType) %>%
    dplyr::group_split() %>%
    purrr::map(function(group_data) {
      sample_split <- dplyr::sample_frac(group_data, ratio_valid)
      list(train = dplyr::anti_join(group_data, sample_split),
           valid = sample_split)
    })

  n <- length(split_data)
  print(n)
  valid_gr <- purrr::map(split_data[1:n], function(split_element) split_element$valid) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(seqnames, start, end) %>%
    GenomicRanges::GRanges()
  train_gr <- purrr::map(split_data[1:n], function(split_element) split_element$train) %>% 
    dplyr::bind_rows() %>%
    dplyr::arrange(seqnames, start, end) %>%
    GenomicRanges::GRanges()

  return(GenomicRanges::GRangesList(train_gr, valid_gr))
}

#' Keep Equal Number of Ranges
keep_equal_num <- function(gr, threshold) {
  random_indices <- sample(length(gr), threshold)
  gr_random <- gr[random_indices]
  gr_random <- GenomicRanges::sort(gr_random)
  return(gr_random)
}

#' Convert thick.start, thick.end, and thick.width to a thick IRanges Column
convert_to_thick_iranges <- function(grl) {
  modified_grl <- lapply(grl, function(gr) {
    # Create IRanges object for the thick column
    thick_ir <- IRanges::IRanges(
      start = mcols(gr)$thick.start,
      end = mcols(gr)$thick.end
    )
    # Assign the IRanges object to the thick column
    S4Vectors::mcols(gr)$thick <- thick_ir
    # Optionally, drop other columns if needed
    S4Vectors::mcols(gr) <- S4Vectors::mcols(gr)[, "thick", drop = FALSE]
    return(gr)
  })

  return(GenomicRanges::GRangesList(modified_grl))
}

#' Wrap Up Data Splitting into Training and Validation Sets
wrapup_split_train_valid <- function(pos_grl, neg_grl) {
  pos_tr_grl <- GenomicRanges::GRangesList()
  pos_vl_grl <- GenomicRanges::GRangesList()

  neg_tr_grl <- GenomicRanges::GRangesList()
  neg_vl_grl <- GenomicRanges::GRangesList()

  for (i in 1:length(pos_grl)) {
    pos_df <- data.frame(pos_grl[[i]])
    neg_df <- data.frame(neg_grl[[i]])

    pos_splt_grl <- split_train_valid(pos_df, ratio_valid = validation_ratio)
    neg_splt_grl <- split_train_valid(neg_df, ratio_valid = validation_ratio)

    final_num_tr <- min(length(pos_splt_grl[[1]]), length(neg_splt_grl[[1]]))
    final_num_vl <- min(length(pos_splt_grl[[2]]), length(neg_splt_grl[[2]]))

    pos_tr_grl[[i]] <- keep_equal_num(pos_splt_grl[[1]], final_num_tr)
    pos_vl_grl[[i]] <- keep_equal_num(pos_splt_grl[[2]], final_num_vl)

    neg_tr_grl[[i]] <- keep_equal_num(neg_splt_grl[[1]], final_num_tr)
    neg_vl_grl[[i]] <- keep_equal_num(neg_splt_grl[[2]], final_num_vl)
  }

  pos_tr_grl <- convert_to_thick_iranges(pos_tr_grl)
  pos_vl_grl <- convert_to_thick_iranges(pos_vl_grl)
  neg_tr_grl <- convert_to_thick_iranges(neg_tr_grl)
  neg_vl_grl <- convert_to_thick_iranges(neg_vl_grl)

  grl_list <- list(
    pos_tr_grl = pos_tr_grl,
    pos_vl_grl = pos_vl_grl,
    neg_tr_grl = neg_tr_grl,
    neg_vl_grl = neg_vl_grl
  )

  return(grl_list)
}

# TESTING: keep equal number
wrapup_test <- function(pos_grl, neg_grl) {
  pos_te_grl <- GenomicRanges::GRangesList()
  neg_te_grl <- GenomicRanges::GRangesList()

  for (i in seq_along(pos_grl)) {
    final_num_te <- min(length(pos_grl[[i]]), length(neg_grl[[i]]))

    pos_te_grl[[i]] <- keep_equal_num(pos_grl[[i]], final_num_te)
    neg_te_grl[[i]] <- keep_equal_num(neg_grl[[i]], final_num_te)
  }

  # make thick column to IRanges
  pos_te_grl <- convert_to_thick_iranges(pos_te_grl)
  neg_te_grl <- convert_to_thick_iranges(neg_te_grl)

  te_ls <- list(pos_te_grl = pos_te_grl,
                neg_te_grl = neg_te_grl)

  return(te_ls)

}



# 1 load Rdata
writeLines("\nLoad data..")
load(ocr_rdata_file)
load(tss_rdata_file)


# 2 select OCR for each library with condition
#In summary, the function filters OCRs based on two main conditions:
# 1. OCRs with more than one-base overlap are always selected.
# 2. Singleton OCRs are selected only if their score is
# greater than or equal to min_score_for_onebase_ocr (default is 4).
writeLines("\nSelect OCR for each library..") # nolint: line_length_linter.
writeLines("\n..")

tr_pos_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocr_train_gr)
tr_pos_aug_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocr_train_aug_gr)

te_pos_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocr_test_gr)
te_pos_aug_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocr_test_aug_gr)

# tr_neg_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocrlike_neg_train_gr) # nolint: line_length_linter.

tr_neg_aug_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, ocrlike_neg_train_aug_gr) # nolint: line_length_linter.
tr_neg_npn_sltocr_grl <- select_ocr_on_condition(tr_orig_ctss, nearpos_neg_train_aug_gr) # nolint: line_length_linter.

te_neg_aug_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocrlike_neg_test_aug_gr) # nolint: line_length_linter.
te_neg_npn_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, nearpos_neg_test_aug_gr) # nolint: line_length_linter.

# te_neg_sltocr_grl <- select_ocr_on_condition(te_orig_ctss, ocrlike_neg_test_gr) # nolint: line_length_linter.



# 3 split data
writeLines("\nSplit data..")

# TRAINING POSITIVE + TESTING POSITIVE
tr_pos_grl <- subset_annotation_positive(tr_pos_sltocr_grl)
tr_pos_aug_grl <- subset_annotation_positive(tr_pos_aug_sltocr_grl)
te_pos_grl <- subset_annotation_positive(te_pos_sltocr_grl)
te_pos_aug_grl <- subset_annotation_positive(te_pos_aug_sltocr_grl)

# TRAINING NEGATIVE + TESTING NEGATIVE
# subsampling negative OCR position to match number of positive,
# but keep annotation ratio of negative
# nearpos negatives were match to positive and
# aug negatives were match to aug positive
tr_neg_npn_grl <- select_negative_to_match_positive(tr_neg_npn_sltocr_grl, tr_pos_grl) # nolint: line_length_linter.
tr_neg_aug_grl <- select_negative_to_match_positive(tr_neg_aug_sltocr_grl, tr_pos_aug_grl) # nolint: line_length_linter.
te_neg_npn_grl <- select_negative_to_match_positive(te_neg_npn_sltocr_grl, te_pos_grl) # nolint: line_length_linter.
te_neg_aug_grl <- select_negative_to_match_positive(te_neg_aug_sltocr_grl, te_pos_aug_grl) # nolint: line_length_linter.

# TRAINING: split training (80%) and validation (20%)
# based on their natural ratio of annotation
tr_vl_ls <- wrapup_split_train_valid(tr_pos_grl, tr_neg_npn_grl)
tr_vl_aug_ls <- wrapup_split_train_valid(tr_pos_aug_grl, tr_neg_aug_grl)


te_ls <- wrapup_test(te_pos_grl, te_neg_npn_grl)
te_aug_ls <- wrapup_test(te_pos_aug_grl, te_neg_aug_grl)


# 4 save
writeLines("\nSaving..")
save(list = c("tr_vl_ls", "tr_vl_aug_ls", "te_ls", "te_aug_ls"),
     file = paste0("5_", name, "_tr_vl_te.RData"))