# Load necessary libraries
suppressPackageStartupMessages({

  library(dotenv)
  library(this.path)

  library(tidyverse)

  library(CAGEfightR)

  library(IRanges)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(GenomicFeatures)

})

# Working directory and environment variables
work_dir <- this.path::this.dir()
setwd(work_dir)
dotenv::load_dot_env("genomewide.env")



## 0. Setting

# Retrieve environment variables
ext_dis <- as.numeric(Sys.getenv("EXT_DIS"))

dir_resources <- Sys.getenv("DIR_RESOURCES")
dir_results <- Sys.getenv("DIR_RESULTS")

infile_ctss_rse <- Sys.getenv("INFILE_CTSS_RSE")
ctss_rse <- readRDS(file.path(dir_resources, infile_ctss_rse))
# if you wanna choose some libraries, you can do it here:
ctss_rse <- ctss_rse[, c(7, 8, 13, 14, 17)] # nolint

outfile_tc_grl <- Sys.getenv("OUTFILE_TC_GRL")

# Genome information and standard chromosomes from UCSC
species <- Sys.getenv("SPECIES")
genome_info <- GenomeInfoDb::keepStandardChromosomes(
  rtracklayer::SeqinfoForUCSCGenome(Sys.getenv("GENOME_INFO")),
  species = species
)
std_chr <- GenomicRanges::seqnames(genome_info)

# Dynamically load TxDb package and create the TxDb object
txdb_package <- Sys.getenv("TXDB_PACKAGE")

if (!requireNamespace(txdb_package, quietly = TRUE)) {
  stop(paste("Package", txdb_package, "not installed"))
}

do.call("library", list(txdb_package))

txdb <- get(txdb_package, envir = asNamespace(txdb_package))



# 1. process the data

ctss_rse <- GenomeInfoDb::keepStandardChromosomes(ctss_rse,
                                                  pruning.mode = "coarse")
col_ctss_rse <- colnames(ctss_rse)

tc_grl <- GRangesList()
for (i in col_ctss_rse){

  ctss <- ctss_rse[, i]
  ctss <- CAGEfightR::calcPooled(ctss, inputAssay = "counts")
  ctss <- subset(ctss, score > 0)

  object <- CAGEfightR::clusterUnidirectionally(ctss)

  new_ranges <- IRanges(start = start(object$thick) - ext_dis,
                        end = end(object$thick) + ext_dis)

  new_object <- object
  ranges(new_object) <- new_ranges

  tc_grl[[i]] <- new_object

}

saveRDS(tc_grl, file.path(dir_results, outfile_tc_grl))
