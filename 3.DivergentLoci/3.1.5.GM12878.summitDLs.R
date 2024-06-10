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
library(data.table)

setwd("/maps/projects/ralab/data/projects/nucleiCAGEproject/3.DivergentLoci/")



divergentLociTCsSummit<- function(ctss, max_gap=400, win_size=200, inputAssay="counts", callingAssay="counts") {
  
  message(paste("Running DL calling for",colnames(ctss)))
  
  ## Split on strand
  ctss.o <- ctss
  ctss <- calcPooled(ctss, inputAssay=callingAssay)
  ctss <- subset(ctss,score>0)
  
  #Call TCs for each sample
  object <- clusterUnidirectionally(ctss)
  object <- quantifyClusters(ctss,object)
  object <- calcPooled(object, inputAssay=inputAssay)
  object <- rowRanges(object) 
  
  message("Removing overlapping TCs by strand...")
  object <- swapRanges(object) #Summit focused
  TCsByStrand <- splitByStrand(object)
  
  ## Find overlapping and book-ended summits between strands
  olaps <- findOverlaps(TCsByStrand$'-',TCsByStrand$'+',maxgap=0,type="any",select="all",ignore.strand=TRUE)
  m_score <-  mcols(TCsByStrand$'-')$score
  p_score <-  mcols(TCsByStrand$'+')$score
  
  m_rem <- queryHits(olaps)[which(m_score[queryHits(olaps)] <= p_score[subjectHits(olaps)])]
  p_rem <- subjectHits(olaps)[which(p_score[subjectHits(olaps)] < m_score[queryHits(olaps)])]
  
  ## remove overlapping TCs
  if (length(m_rem)>0) {
    TCsByStrand$'-' <- TCsByStrand$'-'[-m_rem]
  }
  if (length(p_rem)>0) {
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
  
  mergeLoci <- function(m, p) {
    
    m.i <- 1
    p.i <- 1
    
    while ((p[p.i,start] < m[m.i,start]) || ((p[p.i,start]-m[m.i,start]) > max_gap)) {
      if ((p[p.i,score] < m[m.i,score]) && (p.i < nrow(p)) && ((p[p.i+1,start]-m[m.i,start]) < max_gap)) {
        p.i <- p.i + 1
      } else {
        m.i <- m.i + 1
      }
    }
    floor(m[m.i,start] + ((p[p.i,start] - m[m.i,start]) / 2))
  }
  
  dt <- data.table(strand=as.character(strand(object)),
                   start=start(object),
                   score=object$score,
                   locus=con$membership,
                   stringsAsFactors=FALSE)
  
  dt <- split(dt, dt$strand)
  setorder(dt$'-', -score, start) ## order by score, then by position (prioritize rightmost)
  setorder(dt$'+', -score, -start) ## order by score, then by position (prioritize leftmost)
  
  loci.m <- split(dt$'-', dt$'-'$'locus')
  loci.p <- split(dt$'+', dt$'+'$'locus')
  
  mid <- bplapply(1:length(loci.m), function(i) mergeLoci(loci.m[[i]], loci.p[[i]]))
  mid <- unlist(mid)
  
  ## Extract seqnames for loci
  div_chr <- con$membership[match(1:con$no,con$membership)]
  div_chr <- as.character(sapply(names(div_chr), function(n) strsplit(n,":")[[1]][1]))
  
  #Revert to counts
  ctss <- ctss.o
  ctss <- calcPooled(ctss,inputAssay=inputAssay)
  ctss <- subset(ctss, score>0)
  
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



register(MulticoreParam(workers=30))
#GM12878
#GM12878.CTSSs <- readRDS("../1.CTSSs/repSub.combined.GM12878.CTSSs.rmSingletons.rds")

#Select relevant samples for speedup
#samp <- colnames(GM12878.CTSSs)[c(1:35)]

#GM12878.DLs <- lapply(colnames(GM12878.CTSSs),function(x) divergentLociTCsSummit(GM12878.CTSSs[,x],callingAssay="counts")) #Maybe run later
#GM12878.DLs.rmS <- lapply(samp,function(x) divergentLociTCsSummit(GM12878.CTSSs[,x],callingAssay="counts.noSingletons"))


#saveRDS(GM12878.DLs, "repSub.GM12878.summitDLs.rds")
#saveRDS(GM12878.DLs.rmS, "repSub.GM12878.summitDLs.rmSingletons.rds")


#rm(GM12878.DLs.rmS,GM12878.CTSSs)


GM12878.CTSSs <- readRDS("../1.CTSSs/poolSub.combined.GM12878.CTSSs.rmSingletons.rds")

samp <- colnames(GM12878.CTSSs)[c(1:8,11)]
#GM12878.DLs <- lapply(colnames(GM12878.CTSSs),function(x) divergentLociTCsSummit(GM12878.CTSSs[,x],callingAssay="counts"))
GM12878.DLs.rmS <- lapply(samp,function(n) {
  message("*** ",n," ***")
  CTSSs.by.chr <- split(GM12878.CTSSs[,n], seqnames(GM12878.CTSSs[,n]))
  do.call("c", lapply(names(CTSSs.by.chr), function(chr) {
    message("* ",chr," *")
    divergentLociTCsSummit(CTSSs.by.chr[[chr]],callingAssay="counts.noSingletons")
  }))
})

#saveRDS(GM12878.DLs, "poolSub.GM12878.summitDLs.rds")
saveRDS(GM12878.DLs.rmS, "poolSub.GM12878.summitDLs.rmSingletons.rds")