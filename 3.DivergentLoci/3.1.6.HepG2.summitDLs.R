library(CAGEfightR)

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
library(Matrix)
library(corrplot)
library(data.table)
library(assertthat)


#Function for pooling replicates
poolReplicates <- function(object, replicates, inputAssay="counts") {
  
  ## dimensions
  assert_that(ncol(object)==length(replicates),
              length(unique(replicates))>1)
  
  
  message("Calculating pooled counts...")
  a <- do.call("cbind", lapply(unique(replicates), function(b) 
    Matrix::rowSums(assay(object,inputAssay)[,replicates==b,drop=FALSE])))
  a <- as(a, "sparseMatrix")
  colnames(a) <- unique(replicates)
  
  message("Create a pasted design matrix...")
  c <- do.call("rbind", lapply(unique(replicates), function(b) 
    apply(colData(object)[replicates==b,], 2 , paste , collapse = "," )))
  
  rse <- SummarizedExperiment(assays=SimpleList(counts=a),
                              rowRanges=rowRanges(object), colData=as.data.frame(c))
  
  rse <- calcTotalTags(rse)
  
  rse
}

#Alternative subsampling to deal with dropped rows
subsampleTarget <- function(object, inputAssay = "counts", target) {
  
  assert_that(methods::is(object, "SummarizedExperiment"),
              inputAssay %in% assayNames(object),
              is.numeric(target), target > 0)
  
  a <- assay(object,inputAssay)
  n <- ncol(a)
  nz <- lapply(1:n, function(i) nonzero(a[,i,drop=FALSE])[,1])
  s <- Matrix::colSums(a)
  d <- unlist(lapply(1:n, function(i) {
    if (s[i]<=target)
      a[nz[[i]],i]
    else
      rbinom(length(nz[[i]]),a[nz[[i]],i], target/s[i])
  }))
  keep <- which(sapply(nz,length)>0)
  
  # Create a matrix to store the subsampled counts, ensuring all rows are retained
  new_matrix <- Matrix::sparseMatrix(
    i = unlist(nz),
    j = unlist(lapply(keep, function(i) rep(i, length(nz[[i]])))),
    x = d, 
    dims = dim(a),  # Preserve the original dimensions
    dimnames = list(rownames(a), colnames(a)[keep])
  )
  
  # Replace the original assay with the new matrix
  assay(object, inputAssay) <- new_matrix
  
  # Recalculate the total count of reads for each sample
  object <- calcTotalTags(object)
  
  object
}


rmSingletons <- function(rse, inputAssay = "counts") {
  rmS <- assay(rse, inputAssay)
  # Find the indices where values are equal to 1
  idx <- which(rmS == 1, arr.ind = TRUE)
  # Replace the values at those indices with 0
  rmS[idx] <- 0
  # Assign the modified assay back to the SummarizedExperiment object
  assay(rse, "counts.noSingletons") <- rmS
  rse
}

#Set working directory
setwd("/maps/projects/ralab/data/CAGE/cellLinesCAGE/STAR_map/bw_files/")





#HepG2 design
plus_files <- list.files(path=".", pattern="*HepG2.*.plus.bw")
plus_files <- grep("_2_|_3_", plus_files, value = TRUE)
minus_files <- list.files(path=".", pattern="*HepG2.*.minus.bw")
minus_files <- grep("_2_|_3_", minus_files, value = TRUE)

#Create a simple design matrix
sample_id <- substr(plus_files, 1,10)
type <- substr(plus_files, 1,7)


design <- data.frame(sample=sample_id, type=type, bw_plus=plus_files, bw_minus=minus_files)
row.names(design) <- sample_id



#Format bw_files
plus_files <- design$bw_plus
minus_files <- design$bw_minus
names(plus_files) <- names(minus_files) <- row.names(design)

bw_plus <- BigWigFileList(plus_files)
bw_minus <- BigWigFileList(minus_files)

#Quantify CTSSs 
register(MulticoreParam(workers=20))

CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       design = design)



#Remove non standard chromosomes 
keep <- seqlevels(CTSSs)[1:23]

CTSSs <- keepSeqlevels(CTSSs, keep, pruning.mode = "coarse")

CTSSs <- CTSSs %>% calcTPM() %>% calcPooled()


#Pool replicates
CTSSs.p <- poolReplicates(CTSSs, colData(CTSSs)$type)

CTSSs.p <- calcTPM(CTSSs.p)

CTSSs.p <- calcPooled(CTSSs.p, inputAssay = "counts")


#Subsample to min seq depths of samples
set.seed(1337)

tar <- min(colData(CTSSs.p)$totalTags)

CTSSs.p <- subsampleTarget(CTSSs.p,target=tar)

CTSSs.p <- calcTPM(CTSSs.p)



#Removing singletons 
CTSSs.p <- rmSingletons(CTSSs.p)


#Run divergent loci calling
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



register(MulticoreParam(workers=25))

#Select only relevant samples
samp <- colnames(CTSSs.p)

hepG2.DLs.rmS <- lapply(samp,function(x) divergentLociTCsSummit(CTSSs.p[,x],callingAssay="counts.noSingletons"))

saveRDS(hepG2.DLs.rmS, "/projects/ralab/data/projects/nucleiCAGEproject/3.DivergentLoci/poolSub.hepG2.summitDLs.rmSingletons.rds")
saveRDS(CTSSs.p, "/projects/ralab/data/projects/nucleiCAGEproject/1.CTSSs/poolSub.hepG2.CTSSs.rmSingletons.rds")

