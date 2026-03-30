#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(rtracklayer)
library(CAGEfightR)
library(PRIME)
library(BiocParallel)
library(dplyr)
library(stringr)
})

input_dir <- "../0.External_resources/FANTOM5/pooled_bw"
output_dir <- "PRIME_FANTOM5_facets_TPM_bw"

dir.create(output_dir)

#BPPARAM <- MulticoreParam(workers = 20)
BPPARAM <- SerialParam()

files <- list.files(
  path = input_dir,
  pattern = "\\.bw$",
  full.names = TRUE
)

plus_files  <- files[str_detect(files, "\\.plus\\.bw$")]
minus_files <- files[str_detect(files, "\\.minus\\.bw$")]

get_facet <- function(x) {
  basename(x) %>%
    str_remove("\\.plus\\.bw$") %>%
    str_remove("\\.minus\\.bw$")
}

plus_df <- data.frame(
  facet = get_facet(plus_files),
  plus  = plus_files,
  stringsAsFactors = FALSE
)

minus_df <- data.frame(
  facet = get_facet(minus_files),
  minus = minus_files,
  stringsAsFactors = FALSE
)

facet_df <- inner_join(plus_df, minus_df, by = "facet")

message("Found ", nrow(facet_df), " facets.")

process_facet <- function(i) {

  facet <- facet_df$facet[i]
  plus_file <- facet_df$plus[i]
  minus_file <- facet_df$minus[i]

  message("[", i, "/", nrow(facet_df), "] ", facet)

  bw_plus <- BigWigFileList(plus_file)
  bw_minus <- BigWigFileList(minus_file)

  names(bw_plus)  <- facet
  names(bw_minus) <- facet

  CTSSs <- quantifyCTSSs(
    plusStrand  = bw_plus,
    minusStrand = bw_minus
  )

  stopifnot(!is.null(colnames(CTSSs)))

  CTSSs <- calcTotalTags(CTSSs)
  CTSSs <- calcTPM(CTSSs)

  ## Temporary directory per worker
  tmp_dir <- tempfile(pattern = paste0("bw_", facet, "_"))
  dir.create(tmp_dir)

  PRIME::exportBw(
    object = CTSSs,
    directory = tmp_dir,
    replicates = "all",
    inputAssay = "TPM",
    splitByStrand = TRUE
  )

  bw_files <- list.files(tmp_dir, pattern = "\\.bw$", full.names = TRUE)

  for (f in bw_files) {
    strand <- if (grepl("\\.plus\\.bw$", f)) "plus" else "minus"
    facet_clean <- gsub("\\s+", ".", facet)
    out_name <- paste0(facet_clean, ".", strand, ".bw")
    file.copy(f, file.path(output_dir, out_name), overwrite = TRUE)
  }

  unlink(tmp_dir, recursive = TRUE)

  return(TRUE)
}

res <- bplapply(
  X = seq_len(nrow(facet_df)),
  FUN = process_facet,
  BPPARAM = BPPARAM
)
