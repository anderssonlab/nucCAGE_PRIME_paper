setwd("/projects/ralab/data/projects/nucleiCAGEproject/5.EPCrisprBenchmark")

suppressPackageStartupMessages({
  library(tidyverse)
  library(CAGEfightR)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(tools)
  
  library(rtracklayer)        
  library(tibble)            
  library(readr)              
  library(forcats)           
  library(tidyr)              
  library(GenomicRanges)      
  library(reshape2) 
  
  library(mymlR, lib="/home/zmk214/zmk214/R_packakges/")
  source("2_profile_functions.R")
  
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
})

# 0. setting and imports 

pos.gr <- readRDS("1.EPC.pos.gr.rds")
neg.gr <- readRDS("1.EPC.neg.notsig.gr.rds")

# k562.s.CTSSs <- readRDS("/projects/ralab/data/projects/nucleiCAGEproject/1.CTSSs/subsampled.k562.CTSSs.rds")
# k562.ps.CTSSs <- readRDS("/projects/ralab/data/projects/nucleiCAGEproject/1.CTSSs/pooled.subsampled.k562.CTSSs.rds")
# k562.rs.CTSSs <- readRDS("/projects/ralab/data/projects/nucleiCAGEproject/1.CTSSs/repSub.combined.k562.CTSSs.rds")
k562.ps.c.CTSSs <- readRDS("/projects/ralab/data/projects/nucleiCAGEproject/1.CTSSs/poolSub.combined.k562.CTSSs.rds")

# profile_output_dir <- "EPC_k562.subsampling.CTSSs"
# profile_output_dir <- "EPC_k562.pooled.subsampling.CTSSs"
profile_output_dir <- "EPC_k562.poolSub.combined.CTSSs"

ext <- 250
ATAC_BP_EXT <- ext
len_vec <- ext * 2 + 1


# 1. create the output dir and subfolder 
current_dir <- getwd()
OUTPUT_DIR_NAME <- profile_output_dir
new_folder_path <- file.path(current_dir, OUTPUT_DIR_NAME)
metadata_path <- file.path(new_folder_path , "metadata")
metadata_subtnorm_path <- file.path(new_folder_path , "metadata_subtnorm")
profiles_path <- file.path(new_folder_path , "profiles")
profiles_subtnorm_path <- file.path(new_folder_path , "profiles_subtnorm")

if (!file.exists(new_folder_path)) {
  
  dir.create(new_folder_path)
  dir.create(metadata_path)
  dir.create(metadata_subtnorm_path)
  dir.create(profiles_path)
  dir.create(profiles_subtnorm_path)
  cat("New folder created:", new_folder_path, "\n")
  
  for (i in c(metadata_path, metadata_subtnorm_path, profiles_path, profiles_subtnorm_path)){
    dir.create(file.path(i, "pos"))
    dir.create(file.path(i, "neg"))
  }  
} else {
  cat("Folder already exists:", new_folder_path, "\n")
}



# 3. profile

profile_each_nofilter <- function(cage_rse, atac_gr, set_dir, set_name){
  
  cage_rse <- sortSeqlevels(cage_rse)
  seqlevels(atac_gr) <- seqlevels(cage_rse)
  seqinfo(atac_gr) <- seqinfo(cage_rse)
  
  for (i in 1:ncol(cage_rse)) { # all columns (samples)
    coln <- colnames(cage_rse[,i])
    profile_nofilter(cage_rse, atac_gr, col=coln, setDir=set_dir, setName=set_name)
  }
  
}
profile_each_nofilter(k562.ps.c.CTSSs, pos.gr, "pos", "pos")
profile_each_nofilter(k562.ps.c.CTSSs, neg.gr, "neg", "neg")

