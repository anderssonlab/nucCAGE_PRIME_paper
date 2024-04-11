#CAGEfighR extensions
setwd("/home/mhf817")
sep <- .Platform$file.sep
cfr.extensions.dir <- "/home/mhf817/CAGEfightR_extensions"
source(paste(cfr.extensions.dir, "utils.R", sep = sep))
source(paste(cfr.extensions.dir, "support.R", sep = sep))
source(paste(cfr.extensions.dir, "cumulative.R", sep = sep))
source(paste(cfr.extensions.dir, "decompose.R", sep = sep))
source(paste(cfr.extensions.dir, "enhancers.R", sep = sep))
source(paste(cfr.extensions.dir, "noise.R", sep = sep))
source(paste(cfr.extensions.dir, "coverage.R", sep = sep))
source(paste(cfr.extensions.dir, "complexity.R", sep = sep))
source(paste(cfr.extensions.dir, "subsample.R", sep = sep))
source(paste(cfr.extensions.dir, "dispersion.R", sep = sep))


#Annotations
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Packages for plotting
library(CAGEfightR)
library(epistack)
library(BiocParallel)
library(dplyr)
library(reshape2)
library(limma)
library(ggplot2)
library(parallel)

setwd("/maps/projects/ralab/data/projects/nucleiCAGEproject/3.DivergentLoci/")





#Enhancer calling function


divergentLociTCs2<- function(ctss, max_gap=400, win_size=200, inputAssay="counts") {
  
  
  message(paste("Running DL calling for",colnames(ctss)))
  ## Split on strand
  ctss <- calcPooled(ctss, inputAssay = inputAssay)
  ctss <- subset(ctss,score>0)
  
  #Call TCs for each sample
  object <- clusterUnidirectionally(ctss)
  object <- quantifyClusters(ctss,object)
  object <- calcPooled(object, inputAssay=inputAssay)
  object <- rowRanges(object)
  
  message("Removing overlapping TCs by strand...")
  TCsByStrand <- splitByStrand(object)
  
  ## Find overlaps between TCs on separate strands
  olaps <- findOverlaps(TCsByStrand$'-',TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)
  m_score <-  mcols(TCsByStrand$'-')$score
  p_score <-  mcols(TCsByStrand$'+')$score
  
  m_rem <- queryHits(olaps)[which(m_score[queryHits(olaps)] <= p_score[subjectHits(olaps)])]
  p_rem <- subjectHits(olaps)[which(p_score[subjectHits(olaps)] < m_score[queryHits(olaps)])]
  
  ## remove overlapping TCs
  if(length(m_rem)>0) {
    TCsByStrand$'-' <- TCsByStrand$'-'[-m_rem]
  }
  if(length(p_rem)>0){
    TCsByStrand$'+' <- TCsByStrand$'+'[-p_rem]
  }
  
  
  message("Finding divergent TC pairs...")
  ## Find divergent TC pairs
  m_pad <- flank(TCsByStrand$'-', width=max_gap, start=TRUE, both=FALSE)
  pairs <- findOverlaps(m_pad,TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)
  
  ## Find connected components of TC pair graphs
  edge_list <- cbind(names(TCsByStrand$'-')[queryHits(pairs)],
                     names(TCsByStrand$'+')[subjectHits(pairs)])
  g <- igraph::graph_from_edgelist(edge_list,directed=FALSE)
  con <- igraph::components(g)
  
  ## Keep only relevant TCs
  object <- object[names(con$membership)]
  
  message("Merging into divergent loci...")
  covByStrand <- splitPooled(methods::as(rowRanges(ctss),"GRanges"))
  
  mergeLoci <- function(minus_pos, plus_pos, minus_score, plus_score)
  {
    mid <- NA
    expr <- 0
    
    for (m in 1:length(minus_pos)) {
      for (p in 1:length(plus_pos)) {
        
        if (plus_pos[p] < minus_pos[m])
          next
        
        if (plus_pos[p] - minus_pos[m] > max_gap)
          break
        
        iter_mid <- round(mean(c(minus_pos[m],plus_pos[p])))
        iter_expr <- sum(minus_score[minus_pos<=iter_mid & (minus_pos >= iter_mid-win_size)]) +
          sum(plus_score[plus_pos>=iter_mid & (plus_pos <= iter_mid+win_size)])
        
        if (iter_expr>expr) {
          mid <- iter_mid
          expr <- iter_expr
        }
      }
    }
    
    mid
  }
  
  
  # Set up parallel processing
  bpparam <- MulticoreParam(80)
  

  # Define a function to be executed in parallel
  parallel_function <- function(c) {
    memb=con$membership[con$membership==c]
    d <- object[names(memb)]
    db <- d
    d <- reduce(d, min.gapwidth = max_gap)                  
    vm <- Views(covByStrand$'-', d[strand(d) == "-"])
    vp <- Views(covByStrand$'+', d[strand(d) == "+"])
    sp <- start(unlist(viewRangeMaxs(vp)))
    sp <- sp[length(sp)]  # If multiple summits, take the rightmost
    sm <- end(unlist(viewRangeMaxs(vm)))
    sm <- sm[1]  # If multiple summits, take the leftmost
    mid <- floor(sm + ((sp - sm) / 2))
    if (sm > sp) {  # if summits causes problem, revert to previous solution
      mid <- mergeLoci(end(db[strand(db) == "-"]), start(db[strand(db) == "+"]), db[strand(db) == "-"]$score, db[strand(db) == "+"]$score)
    }
    return(mid)
  }
  

  # Parallelize the function execution
  mid <- bplapply(unique(con$membership),parallel_function, BPPARAM = bpparam)
  
  mid <- unlist(mid)
  ## Extract seqnames for loci
  div_chr <- con$membership[match(1:con$no,con$membership)]
  div_chr <- as.character(sapply(names(div_chr), function(n) strsplit(n,":")[[1]][1]))
  
  gr <- GRanges(seqnames=div_chr,IRanges(start=mid,end=mid))
  seqlevels(ctss,pruning.mode="coarse") <- seqlevels(gr)
  seqinfo(ctss) <- seqinfo(gr)
  
  cat("\r")
  message("Calculating directionality...")
  
  win_1 <- flank(gr,width=win_size,start=TRUE,both=FALSE)
  win_2 <- flank(gr,width=win_size,start=FALSE,both=FALSE)
  
  covByStrand <- splitPooled(methods::as(rowRanges(ctss),"GRanges"))
  ## Quantify strand-wise in flanking windows around midpoint
  M1 <- unlist(viewSums(Views(covByStrand$`-`, win_1)))
  P2 <- unlist(viewSums(Views(covByStrand$`+`, win_2)))
  M2 <- unlist(viewSums(Views(covByStrand$`-`, win_2)))
  P1 <- unlist(viewSums(Views(covByStrand$`+`, win_1)))
  
  ## Calculate directionality
  pooled_directionality <- (P2-M1) / (P2+M1)
  
  ## Test if divergent
  divergent <- (M1>P1) & (P2>M2)
  
  message("Calculating coverage across samples...")
  
  ## Quantify strand-wise in flanking windows around midpoint
  strand(win_1) <- "-"
  strand(win_2) <- "+"
  mat_2_plus <- suppressMessages(assay(quantifyClusters(ctss, win_2, inputAssay = inputAssay),inputAssay) > 0)
  mat_1_minus <- suppressMessages(assay(quantifyClusters(ctss, win_1, inputAssay = inputAssay),inputAssay) > 0)
  
  ## Quntify number of bidirectional cases (both strands expressed)
  bidirectional <- rowSums(mat_1_minus & mat_2_plus)
  
  message("Preparing output...")
  
  ## Build GRanges object
  start(gr) <- start(gr)-win_size
  end(gr) <- end(gr)+win_size
  gr$score <- M1+P2
  gr$thick <- IRanges(start=mid,width=1)
  
  mcols(gr)[, "directionality"] <- pooled_directionality
  mcols(gr)[, "bidirectionality"] <- bidirectional
  mcols(gr)[, "divergent"] <- divergent
  
  ids <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))
  names(gr) <- ids
  names(gr$thick) <- ids
  
  gr$sample <- colnames(ctss)
  ## Remove non-divergent cases
  gr
}


#k562.s.CTSSs <- readRDS("../1.CTSSs/repSub.combined.k562.CTSSs.rds")
#k562.DLs.CTSS <- lapply(colnames(k562.s.CTSSs),function(x) divergentLociTCs2(k562.s.CTSSs[,x]))

#saveRDS(k562.DLs.CTSS,"k562.summitDLs.rds")


k562.s.CTSSs <- readRDS("../1.CTSSs/poolSub.combined.k562.CTSSs.rds")
k562.DLs.CTSS <- lapply(colnames(k562.s.CTSSs),function(x) divergentLociTCs2(k562.s.CTSSs[,x]))

saveRDS(k562.DLs.CTSS,"k562.PoolsummitDLs.rds")

