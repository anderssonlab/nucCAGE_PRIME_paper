#!/usr/bin/env Rscript

writeLines("\n #### 6 CREATE PROFILES ####")

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

### ARGPARSE
parser <- ArgumentParser()
parser$add_argument("-r", "--ocr_rdata", default = ".",
                    help = "OCR and OCR-like RData file ")
parser$add_argument("-c", "--cage_tss_rdata", default = ".",
                    help = "CAGE TSS rData file")
parser$add_argument("-o", "--output_dir", default = ".",
                    help = "output directory")
parser$add_argument("-O", "--output_dir_name", default = "profiles_result",
                    help = "output dir name")
parser$add_argument("-d", "--distance", default = 200, type = "integer",
                    help = "distance from the dhs peak (bps)")
parser$add_argument("-s", "--save_count_profiles", action = "store_true",
                    default = FALSE,
                    help = "Flag to save count profile. Default is FALSE.")

args <- parser$parse_args()

# input and output
ocr_rdata_file <- args$ocr_rdata
tss_rdata_file <- args$cage_tss_rdata
outdir_dir <- args$output_dir
outdir_dir_name <- args$output_dir_name

# Fixed
ext_dis <- as.numeric(args$distance)
save_count_profiles <- args$save_count_profiles

outdir_main_name <- c("metadata",
                      "profiles",
                      "profiles_subtnorm",
                      "predictions")
outdir_subdir_name <- c("pos", "neg")



# 0 helper functions
source("functions.r")

#' Prepare Directory Structure for Profile Output
prep_profile_dir <- function(output_dir = ".",
                             output_dir_name = "profile_output",
                             output_main_name = c("metadata",
                                                  "profiles",
                                                  "profiles_subtnorm",
                                                  "predictions"),
                             output_subdir_name = c("pos", "neg")) {

  # Create the output dir and main dir
  new_path <- file.path(output_dir, output_dir_name)

  if (!file.exists(new_path)) {
    dir.create(new_path)

    # Create main directories and their subdirectories
    lapply(output_main_name, function(main) {
      main_path <- file.path(new_path, main)
      dir.create(main_path)

      if (length(output_subdir_name) > 0) {
        lapply(output_subdir_name, function(subdir) {
          subdir_path <- file.path(main_path, subdir)
          dir.create(subdir_path)
        })
      }
    })

    cat("New folder created:", new_path, "\n")
  } else {
    cat("Folder already exists:", new_path, "\n")
  }
}

#' Normalized Strand Subtraction
strands_norm_subtraction <- function(vec, len_vec) {
  p <- as.numeric(vec[1:len_vec])
  m <- as.numeric(vec[(len_vec + 1):(len_vec * 2)])
  # Normalized strand subtraction
  return((p - m) / max(abs(c(p, m))))
}

#' Apply Normalized Strand Subtraction to All Windows
strands_norm_subtraction_all <- function(windows, ext_dis, len_vec) {
  #' Apply normalized forward and reverse subtraction to all windows
  #'
  p_min_m_norm_df <- as.data.frame(t(apply(windows, 1,
                                           strands_norm_subtraction, len_vec)))
  pos <- seq(1, ncol(p_min_m_norm_df))
  colnames(p_min_m_norm_df) <- paste0("Pos", pos - ext_dis - 1)
  return(p_min_m_norm_df)
}

#' Wrap up profile for input of the model
#'
#' @param ctss_rse SummarizedExperiment object containing count data.
#' @param regions_gr GRanges or GRangesList object specifying genomic regions.
#' @param output_dir Directory where output profiles and metadata will be saved.
#' @param output_dir_name Output directory name within \code{dir_results}.
#' @param output_subdir_name Subdirectory name within \code{outdir_dir_name}
#' where output files will be stored.
#' @param ext_dis Numeric value specifying the extension distance.
#' @param addtn_to_filename Additional text to append to the output file names.
#' @param save_count_profiles Logical value indicating
#' whether to save count profiles.
#' @param file_type Character string indicating the file type for output 
#' ("csv" or "parquet"). Default is "parquet".
wrapup_make_profiles <- function(ctss_rse,
                                 regions_gr,
                                 output_dir,
                                 output_dir_name,
                                 output_subdir_name,
                                 ext_dis,
                                 addtn_to_filename = "",
                                 save_count_profiles = FALSE,
                                 file_type = "parquet") {
  for (i in seq_along(SummarizedExperiment::colnames(ctss_rse))) {

    print(SummarizedExperiment::colnames(ctss_rse)[i])

    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
    print(formatted_datetime)

    if (inherits(regions_gr, "GRangesList")) {
      current_region_gr <- regions_gr[[i]]
    } else if (inherits(regions_gr, "GRanges")) {
      current_region_gr <- regions_gr
    } else {
      stop("regions_gr is neither GRanges nor GRangesList")
    }

    ctss_gr <- cast_rse_to_granges(ctss_rse, assay = "counts", coln_assay = i)

    current_region_gr <- convert_strand_to_nostrand_gr(current_region_gr)
    current_region_gr <- remove_metadata_and_duplicates(current_region_gr)

    print("Start making profile")
    count_profiles <- PRIME::heatmapData(current_region_gr, ctss_gr)
    print("Finish making profile")
    rm(current_region_gr, ctss_gr)

    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
    print(formatted_datetime)

    len_vec <- ext_dis * 2 + 1

    combined_count_profiles <- combine_plus_minus_profiles(count_profiles, len_vec)
    rm(count_profiles)

    combined_subtnorm_profiles <- strands_norm_subtraction_all(combined_count_profiles, ext_dis, len_vec)

    combined_count_metadata <- create_granges_from_rownames(rownames(combined_count_profiles))
    sum_count <- data.frame(rowSums(combined_count_profiles))
    colnames(sum_count) <- "sum_count"
    combined_count_metadata$sum_count <- sum_count

    # Add rownames as a column before saving
    combined_count_metadata$rownames <- rownames(combined_count_metadata)

    if (file_type == "csv") {
      # Save as CSV without rownames
      write.csv(as.data.frame(combined_count_metadata),
                file = file.path(output_dir,
                                 output_dir_name,
                                 "metadata",
                                 output_subdir_name,
                                 paste0("metadata_count_", output_subdir_name,
                                        addtn_to_filename, "_",
                                        SummarizedExperiment::colnames(ctss_rse)[i], ".csv")),
                row.names = FALSE)
    } else if (file_type == "parquet") {
      # Save as Parquet
      arrow::write_parquet(as.data.frame(combined_count_metadata),
                           file.path(output_dir,
                                     output_dir_name,
                                     "metadata",
                                     output_subdir_name,
                                     paste0("metadata_count_", output_subdir_name,
                                            addtn_to_filename, "_",
                                            SummarizedExperiment::colnames(ctss_rse)[i], ".parquet")))
    }
    rm(combined_count_metadata)

    if (save_count_profiles) {
      combined_count_profiles$rownames <- rownames(combined_count_profiles)

      if (file_type == "csv") {
        # Save as CSV without rownames
        write.csv(as.data.frame(combined_count_profiles),
                  file = file.path(output_dir,
                                   output_dir_name,
                                   "profiles",
                                   output_subdir_name,
                                   paste0("profiles_count_", output_subdir_name,
                                          addtn_to_filename, "_",
                                          SummarizedExperiment::colnames(ctss_rse)[i], ".csv")),
                  row.names = FALSE)
      } else if (file_type == "parquet") {
        # Save as Parquet
        arrow::write_parquet(as.data.frame(combined_count_profiles),
                             file.path(output_dir,
                                       output_dir_name,
                                       "profiles",
                                       output_subdir_name,
                                       paste0("profiles_count_", output_subdir_name,
                                              addtn_to_filename, "_",
                                              SummarizedExperiment::colnames(ctss_rse)[i], ".parquet")))
      }
    }
    rm(combined_count_profiles)

    combined_subtnorm_profiles$rownames <- rownames(combined_subtnorm_profiles)

    if (file_type == "csv") {
      # Save as CSV without rownames
      write.csv(as.data.frame(combined_subtnorm_profiles),
                file = file.path(output_dir,
                                 output_dir_name,
                                 "profiles_subtnorm",
                                 output_subdir_name,
                                 paste0("profiles_subtnorm_", output_subdir_name,
                                        addtn_to_filename, "_",
                                        SummarizedExperiment::colnames(ctss_rse)[i], ".csv")),
                row.names = FALSE)
    } else if (file_type == "parquet") {
      # Save as Parquet
      arrow::write_parquet(as.data.frame(combined_subtnorm_profiles),
                           file.path(output_dir,
                                     output_dir_name,
                                     "profiles_subtnorm",
                                     output_subdir_name,
                                     paste0("profiles_subtnorm_", output_subdir_name,
                                            addtn_to_filename, "_",
                                            SummarizedExperiment::colnames(ctss_rse)[i], ".parquet")))
    }
    rm(combined_subtnorm_profiles)

    gc()

    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
    print(current_datetime)
  }
}



# 1 load Rdata
load(ocr_rdata_file)
load(tss_rdata_file)

# 2 create the output dir
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)

# 3 create profile and its metadata

# Define a function to create profiles for a single subdir_name
report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$pos_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_tr",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$pos_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_vl",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

writeLines("DONE POSITIVE..")

report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$neg_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_tr",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss, tr_vl_ls$neg_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_vl",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

writeLines("DONE NEGATIVE..")

report_time_execution(wrapup_make_profiles(te_orig_ctss, te_ls$pos_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_te",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(te_orig_ctss, te_ls$neg_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_te",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

report_time_execution(wrapup_make_profiles(te_orig_ctss, te_aug_ls$pos_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_te_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(te_orig_ctss, te_aug_ls$neg_te_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_te_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))

writeLines("DONE TESTING..")

report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$pos_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_tr_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$pos_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "pos", ext_dis,
                                           addtn_to_filename = "_vl_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$neg_tr_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_tr_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
report_time_execution(wrapup_make_profiles(tr_orig_ctss,
                                           tr_vl_aug_ls$neg_vl_grl,
                                           outdir_dir, outdir_dir_name,
                                           "neg", ext_dis,
                                           addtn_to_filename = "_vl_aug",
                                           save_count_profiles = TRUE,
                                           file_type = "csv"))
writeLines("DONE AUGMENTATION..")