setwd("/maps/projects/ralab/people/mhf817/PL_dcTapSeq")

library(CAGEfightR)
library(GenomicRanges)
library(readr)
library(rtracklayer)
library(parallel)
library(ggplot2)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SummarizedExperiment)
library(tibble)
library(ggrepel)
library(tidyr)

pooled.PRIMEloci.files <- c("K562_sld_C" = "/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/K562_sld_on_PRIMEloci1.0_C.bed",
                            "K562_sld_N" = "/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/K562_sld_on_PRIMEloci1.0_N.bed",
                            "HepG2_sld_C" = "/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/HepG2_sld_on_PRIMEloci1.0_C.bed",
                            "HepG2_sld_N" = "/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/HepG2_sld_on_PRIMEloci1.0_N.bed",
                            "GM12878_sld_C" ="/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/GM12878_sld_on_PRIMEloci1.0_C_1M.bed",
                            "GM12878_sld_N" = "/projects/ralab/data/projects/nucleiCAGEproject/8.Genomewide_prediction/GM12878_sld_on_PRIMEloci1.0_N_5M.bed")


pooled.PRIMEloci <- lapply(names(pooled.PRIMEloci.files), function(n) {
  gr <- read_tsv(pooled.PRIMEloci.files[n],show_col_types = FALSE)
  score <- gr$score
  gr <- makeGRangesFromDataFrame(gr[,1:3])
  gr$score <- score
  gr <- resize(gr, width=200, fix="center")
  gr
})

names(pooled.PRIMEloci) <- names(pooled.PRIMEloci.files)

dc_tap <- read_tsv("media-3.tsv")
dc_tap_hg38 <- dc_tap[,2:4]
colnames(dc_tap_hg38) <- c("chr","start","end")
dc_tap_hg38 <- makeGRangesFromDataFrame(dc_tap_hg38)


# Find overlaps
olap <- findOverlaps(dc_tap_hg38, pooled.PRIMEloci[["K562_sld_N"]])

# If no overlaps at all, just add NA columns and return
if (length(olap) == 0L) {
  mcols(dc_tap_hg38)$score_max_k562     <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_min_k562     <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_mean_k562    <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_median_k562  <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$midpoint_k562      <- rep(NA_real_, length(dc_tap_hg38))
} else {
  subj <- pooled.PRIMEloci[["K562_sld_N"]]
  
  # group scores by query index
  scores_by_query <- split(mcols(subj)$score[subjectHits(olap)], queryHits(olap))
  
  # compute midpoints of overlapping ranges
  midpoints <- start(subj) + (width(subj) %/% 2)
  mids_by_query <- split(midpoints[subjectHits(olap)], queryHits(olap))
  
  n <- length(dc_tap_hg38)
  score_max     <- rep(NA_real_, n)
  score_min     <- rep(NA_real_, n)
  score_mean    <- rep(NA_real_, n)
  score_median  <- rep(NA_real_, n)
  midpoint      <- rep(NA_real_, n)
  
  # safe summarizer: returns NA if no usable values
  safe_apply <- function(x, FUN) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0L) return(NA_real_)
    FUN(x2)
  }
  
  # compute summaries
  vals_max     <- sapply(scores_by_query, function(x) safe_apply(x, max))
  vals_min     <- sapply(scores_by_query, function(x) safe_apply(x, min))
  vals_mean    <- sapply(scores_by_query, function(x) safe_apply(x, mean))
  vals_median  <- sapply(scores_by_query, function(x) safe_apply(x, median))
  
  # for midpoint, usually one value per overlap; if multiple, take mean
  vals_mid     <- sapply(mids_by_query, function(x) safe_apply(x, mean))
  
  # fill in
  idx <- as.integer(names(scores_by_query))
  score_max[idx]     <- vals_max
  score_min[idx]     <- vals_min
  score_mean[idx]    <- vals_mean
  score_median[idx]  <- vals_median
  midpoint[idx]      <- vals_mid
  
  # attach
  mcols(dc_tap_hg38)$score_max_k562     <- score_max
  mcols(dc_tap_hg38)$score_min_k562     <- score_min
  mcols(dc_tap_hg38)$score_mean_k562    <- score_mean
  mcols(dc_tap_hg38)$score_median_k562  <- score_median
  mcols(dc_tap_hg38)$midpoint_k562      <- midpoint
}


# Find overlaps
olap <- findOverlaps(dc_tap_hg38, pooled.PRIMEloci[["HepG2_sld_N"]])


# If no overlaps at all, just add NA columns and return
if (length(olap) == 0L) {
  mcols(dc_tap_hg38)$score_max    <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_min    <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_mean   <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_median <- rep(NA_real_, length(dc_tap_hg38))
} else {
  # group scores of pooled.PRIMEloci[["PRIMEloci CAGE cells pooled"]] by query index of dc_tap_hg38
  scores_by_query <- split(mcols(pooled.PRIMEloci[["HepG2_sld_N"]])$score[subjectHits(olap)], queryHits(olap))
  
  n <- length(dc_tap_hg38)
  score_max    <- rep(NA_real_, n)
  score_min    <- rep(NA_real_, n)
  score_mean   <- rep(NA_real_, n)
  score_median <- rep(NA_real_, n)
  
  # safe summarizer: returns NA if all values are NA / length 0
  safe_apply <- function(x, FUN) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0L) return(NA_real_)
    FUN(x2)
  }
  
  # compute summaries for only the queries that had overlaps
  vals_max    <- sapply(scores_by_query, function(x) safe_apply(x, max))
  vals_min    <- sapply(scores_by_query, function(x) safe_apply(x, min))
  vals_mean   <- sapply(scores_by_query, function(x) safe_apply(x, mean))
  vals_median <- sapply(scores_by_query, function(x) safe_apply(x, median))
  
  # fill into the full-length vectors using the query indices (names of the list)
  idx <- as.integer(names(scores_by_query))
  score_max[idx]    <- vals_max
  score_min[idx]    <- vals_min
  score_mean[idx]   <- vals_mean
  score_median[idx] <- vals_median
  
  # attach to dc_tap_hg38
  mcols(dc_tap_hg38)$score_max_hepg2    <- score_max
  mcols(dc_tap_hg38)$score_min_hepg2    <- score_min
  mcols(dc_tap_hg38)$score_mean_hepg2   <- score_mean
  mcols(dc_tap_hg38)$score_median_hepg2 <- score_median
}


# Find overlaps
olap <- findOverlaps(dc_tap_hg38, pooled.PRIMEloci[["GM12878_sld_N"]])


# If no overlaps at all, just add NA columns and return
if (length(olap) == 0L) {
  mcols(dc_tap_hg38)$score_max    <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_min    <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_mean   <- rep(NA_real_, length(dc_tap_hg38))
  mcols(dc_tap_hg38)$score_median <- rep(NA_real_, length(dc_tap_hg38))
} else {
  # group scores of pooled.PRIMEloci[["PRIMEloci CAGE cells pooled"]] by query index of dc_tap_hg38
  scores_by_query <- split(mcols(pooled.PRIMEloci[["GM12878_sld_N"]])$score[subjectHits(olap)], queryHits(olap))
  
  n <- length(dc_tap_hg38)
  score_max    <- rep(NA_real_, n)
  score_min    <- rep(NA_real_, n)
  score_mean   <- rep(NA_real_, n)
  score_median <- rep(NA_real_, n)
  
  # safe summarizer: returns NA if all values are NA / length 0
  safe_apply <- function(x, FUN) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0L) return(NA_real_)
    FUN(x2)
  }
  
  # compute summaries for only the queries that had overlaps
  vals_max    <- sapply(scores_by_query, function(x) safe_apply(x, max))
  vals_min    <- sapply(scores_by_query, function(x) safe_apply(x, min))
  vals_mean   <- sapply(scores_by_query, function(x) safe_apply(x, mean))
  vals_median <- sapply(scores_by_query, function(x) safe_apply(x, median))
  
  # fill into the full-length vectors using the query indices (names of the list)
  idx <- as.integer(names(scores_by_query))
  score_max[idx]    <- vals_max
  score_min[idx]    <- vals_min
  score_mean[idx]   <- vals_mean
  score_median[idx] <- vals_median
  
  # attach to dc_tap_hg38
  mcols(dc_tap_hg38)$score_max_gm12878    <- score_max
  mcols(dc_tap_hg38)$score_min_gm12878    <- score_min
  mcols(dc_tap_hg38)$score_mean_gm12878   <- score_mean
  mcols(dc_tap_hg38)$score_median_gm12878 <- score_median
}

dc_tap_hg38_df <- cbind(as.data.frame(dc_tap_hg38), dc_tap[,9:ncol(dc_tap)])

ggplot(dc_tap_hg38_df, aes(dist<=15, log10(abs(pct_change_effect_size)),fill=cell_type)) + geom_boxplot() 

ggplot(dc_tap_hg38_df, aes(score_max, fold_change_effect_size,color=cell_type)) +geom_point() + theme_bw()

ggplot(dc_tap_hg38_df, aes(score_max, log2(fold_change_effect_size),color=element_category)) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(score_max, log2(fold_change_effect_size),color=element_location)) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(score_max, log2(fold_change_effect_size),color=design_file_type)) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(score_max, -log10(sceptre_adj_p_value),color=element_location)) +geom_point() + theme_bw()


ggplot(dc_tap_hg38_df, aes(score_max, log2(fold_change_effect_size),color=cell_type)) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(element_category, score_max)) + geom_boxplot() + theme_bw() 
ggplot(dc_tap_hg38_df, aes(significant, score_max,fill=element_location)) + geom_boxplot() + theme_bw()+ facet_wrap(~cell_type)
ggplot(dc_tap_hg38_df, aes(significant, score_max,fill=cell_type)) + geom_boxplot() + theme_bw()


ggplot(dc_tap_hg38_df, aes(score_max_k562, pct_change_effect_size,color=log10(1+dist))) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(score_max_hepg2, pct_change_effect_size,color=significant)) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(score_max_gm12878, pct_change_effect_size,color=significant)) +geom_point() + theme_bw()
ggplot(dc_tap_hg38_df, aes(dist, pct_change_effect_size)) +geom_point() + theme_bw()+ facet_wrap(~element_category)


ggplot(dc_tap_hg38_df, aes(dist, score_max_k562,color=log10(abs(pct_change_effect_size)))) +geom_point() + theme_bw() + facet_wrap(~element_category) + scale_color_gradient2()


