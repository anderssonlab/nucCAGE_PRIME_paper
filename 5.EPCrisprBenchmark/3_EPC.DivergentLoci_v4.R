setwd("/projects/ralab/data/projects/nucleiCAGEproject/5.EPCrisprBenchmark")

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
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
  
  library(argparse)
  require(GenomicFeatures)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  library(tools)
  require(Matrix)
  require(caTools)
  require(assertthat)
  require(ggplot2)
})  


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#### DL

# 1. prep DL, select only bidirectional loci and enhancers
#subsampled.k562.DLs.TCs <- readRDS("/projects/ralab/data/projects/nucleiCAGEproject/3.DivergentLoci/subsampled.k562.DLs.TCs.rds")
#subsampled.k562.DLs.TCs <- keepStandardChromosomes(subsampled.k562.DLs.TCs, pruning.mode = "coarse")

repSub.k562.DLs.TCs <- readRDS("/projects/ralab/data/projects/nucleiCAGEproject/3.DivergentLoci/repSub.k562.DLs.TCs.rds")
repSub.k562.DLs.TCs <- keepStandardChromosomes(repSub.k562.DLs.TCs, pruning.mode = "coarse")

## Filter for rows for positive DLs
DL.gr <- repSub.k562.DLs.TCs[repSub.k562.DLs.TCs$divergent == TRUE &
                               repSub.k562.DLs.TCs$bidirectionality == 1 &
                               repSub.k562.DLs.TCs$score >= 3, ]
# &repSub.k562.DLs.TCs$directionality < 0.8

## Build negative filter (mask) for enhancers
exons <- exons(txdb)
mcols(exons) <- NULL
promoters <- promoters(txdb)
mcols(promoters) <- NULL
neg.filter <- reduce(c(exons,promoters))
## enhancer 
DL.enh.gr <- DL.gr[!(DL.gr %over% neg.filter)]



# 2. import EPC object
load("1.EPC.DHS.eQTL.RData")



# 3. find overlap between EPC and DLs/enh

EPCoverlapDLs <- function(epc.gr, dl.gr, sample.name){
  # find overlap between EPC and chosen DLs object
  ovl <- findOverlaps(epc.gr, dl.gr[(dl.gr$sample == sample.name),]) 
  
  # get directionality of overlaped DLs
  directionality <- c()
  for(i in 1:length(ovl)){
    tmp <- dl.gr[subjectHits(ovl)[i]]$directionality %>% as.numeric(unlist(.)) 
    directionality <- c(directionality, tmp)
  }
  
  # if the same EPC row got multiple hits for DLs,
  # select the directionality that has LOWEST absolute directionality
  # some EPC got same abs directionality, but diff direction (+, -)
  # that's why I keep the abs directionality
  # distinct = no duplicate
  # keep abs directionality for further analysis
  
  ovl.df <- ovl %>% data.frame()
  ovl.df$directionality <- directionality
  ovl.df$abs.directionality <- abs(directionality)
  slt.ovl.df <- ovl.df %>% dplyr::select(queryHits, abs.directionality)
  slt.ovl.df <- slt.ovl.df %>% dplyr::group_by(queryHits) %>% 
    dplyr::filter(abs.directionality == min(abs.directionality)) %>% 
    dplyr::distinct(., queryHits, abs.directionality) %>%
    data.frame()
  
  # if no overlap between EPC and DLs, abs directionality keep NA 
  # sort by queryHits(EPC)
  vec <- c(1:length(epc.gr))
  vec1 <- vec[!vec %in% slt.ovl.df$queryHits] 
  vec2 <- rep(NA, length(vec1))
  vec.df <- data.frame(queryHits=vec1, 
                       abs.directionality=vec2,
                       check.names = FALSE)        # Prevent automatic adjustment of column names
  
  abs.dir <- rbind(slt.ovl.df, vec.df) 
  abs.dir <- abs.dir[order(abs.dir$queryHits), ]
  
  return(abs.dir$abs.directionality)
}

EPCdirectionality <- function(epc.gr, dl.gr){
  ls.samples <- dl.gr$sample %>% unique()
  all.dir.df = c()
  for(i in ls.samples){
    print(i)
    abs.directionality <- EPCoverlapDLs(epc.gr, dl.gr, i)
    print(head(abs.directionality))
    coln = paste0('absdir.', i)
    all.dir.df[[coln]] <- abs.directionality
  }
  all.dir.df <- data.frame(all.dir.df)
  return(all.dir.df)
}  

wrapup_EPCdirectionality <- function(epc.gr, dl.gr, filename, threshold){
  df <- EPCdirectionality(epc.gr, dl.gr) %>% data.frame()
  rownames(df) <- names(epc.gr)
  # the csv here will be the abs.dir score
  write.csv(df, filename)
  
  # comment this line if wanna get the abs.dir value instead
  df <- data.frame(sapply(df, function(x) ifelse(is.na(x), NA, x < threshold))) 
  epc.df <- epc.gr %>% data.frame()
  epc.DL.df <- cbind(epc.df, df)
  
  epc.DL.gr <- epc.DL.df %>% GRanges()
  
  return(epc.DL.gr)
}

pos.DL.gr <- wrapup_EPCdirectionality(pos.gr, DL.gr, "3.pos.DL.csv", 0.8)
# pos.DL.enh.gr <- wrapup_EPCdirectionality(pos.gr, DL.enh.gr, "3.pos.DL.enh.csv", 0.8)

neg.DL.gr <- wrapup_EPCdirectionality(neg.notsig.gr, DL.gr, "3.neg.DL.csv", 0.8)
# neg.DL.enh.gr <- wrapup_EPCdirectionality(neg.notsig.gr, DL.enh.gr, "3.neg.DL.enh.csv", 0.8)

save(list=c("pos.DL.gr", "neg.DL.gr"), file=paste0("3.EPCxDL.DHS.eQTL.RData"))
