#!/usr/bin/env Rscript
# cellline_bigwigs_TPM.R
#
# Generate pooled TPM-normalized bigWig tracks for cell lines
# (GM12878, K562, HCT116, HepG2, A549) using saved CTSS RDS objects
# from 1.CTSSs/ and PRIME::writeBw.
#
# Pools are defined by the Type column in each design matrix.
# Output: PRIME_cellLines_TPM_bw/

suppressPackageStartupMessages({
    library(CAGEfightR)
    library(PRIME)
    library(SummarizedExperiment)
    library(BiocGenerics)
})

ctss_dir   <- "../1.CTSSs"
output_dir <- "PRIME_cellLines_TPM_bw"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: ensure TPM assay is present ----
ensure_tpm <- function(obj) {
    if (!"TPM" %in% assayNames(obj)) {
        obj <- CAGEfightR::calcTotalTags(obj)
        obj <- CAGEfightR::calcTPM(obj)
    }
    obj
}

# ---- Helper: write bigWigs for all samples in a CTSS object ----
write_bw_samples <- function(obj, output_dir) {
    obj <- ensure_tpm(obj)
    for (nm in colnames(obj)) {
        message("  Writing bigWig: ", nm)
        PRIME::writeBw(
            object        = obj,
            dir           = output_dir,
            replicates    = nm,
            inputAssay    = "TPM",
            splitByStrand = TRUE
        )
    }
}

# ---- Helper: load design matrix with flexible column names ----
load_design <- function(cell_line) {
    design_file <- file.path(ctss_dir, paste0(cell_line, ".design.txt"))
    design <- read.table(design_file, sep = "\t", header = TRUE,
                         check.names = FALSE)
    # Filter to CAGE method when column is present
    if ("method" %in% colnames(design)) {
        design <- design[design$method == "CAGE", ]
    }
    # Filter to use==yes when column is present
    if ("use" %in% colnames(design)) {
        design <- design[design$use == "yes", ]
    }
    design
}

# ============================================================
# K562
# ============================================================
message("Processing K562...")
{
    design <- load_design("K562")
    cage_types <- unique(design$Type)

    rds_file <- file.path(ctss_dir, "all.combined.k562.CTSSs.rds")
    ctss <- readRDS(rds_file)

    # The combined object contains individual replicates (pooled=="no") and
    # pooled samples (pooled=="yes"). Use the pooled layer which has one
    # column per (Cell.line, Type) combination.
    cd <- colData(ctss)
    # Pooled colnames are "K562_<Type>"; keep only CAGE types.
    cage_colnames <- paste0("K562_", cage_types)
    keep <- colnames(ctss) %in% cage_colnames
    ctss_cage <- ctss[, keep]

    write_bw_samples(ctss_cage, output_dir)
}

# ============================================================
# GM12878
# ============================================================
message("Processing GM12878...")
{
    rds_file <- file.path(ctss_dir, "all.combined.GM12878.CTSSs.rds")
    ctss <- readRDS(rds_file)

    # The combined object has:
    #   pooled=="yes"  -> pooled by (Cell.line, cell_count, Type)
    #   pooled=="no"   -> pooled by ID (across sequencing runs)
    # We want the cell_count+Type pools (pooled=="yes") that are flagged
    # as part of the main analysis (main=="yes"), which excludes low-input
    # amounts (50K, 250K), triptolide, and GRO-cap samples.
    cd <- colData(ctss)
    keep <- cd$pooled == "yes" & cd$main == "yes"
    ctss_main <- ctss[, keep]

    write_bw_samples(ctss_main, output_dir)
}

# ============================================================
# HCT116
# ============================================================
message("Processing HCT116...")
{
    rds_file <- file.path(ctss_dir, "poolSub.HCT116.CTSSs.rmSingletons.rds")
    ctss <- readRDS(rds_file)

    # Prefix colnames with cell line name
    colnames(ctss) <- paste0("HCT116_", colnames(ctss))

    write_bw_samples(ctss, output_dir)
}

# ============================================================
# HepG2
# ============================================================
message("Processing HepG2...")
{
    rds_file <- file.path(ctss_dir, "poolSub.HepG2.CTSSs.rmSingletons.rds")
    ctss <- readRDS(rds_file)

    colnames(ctss) <- paste0("HepG2_", colnames(ctss))

    write_bw_samples(ctss, output_dir)
}

# ============================================================
# A549
# ============================================================
message("Processing A549...")
{
    rds_file <- file.path(ctss_dir, "poolSub.A549.CTSSs.rmSingletons.rds")
    ctss <- readRDS(rds_file)

    colnames(ctss) <- paste0("A549_", colnames(ctss))

    write_bw_samples(ctss, output_dir)
}

message("Done. BigWigs written to: ", output_dir)
