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





##################################

#Summit based enhancer calling


divergentLociTCsSummit<- function(ctss, max_gap=400, win_size=200, inputAssay="counts", callingAssay="counts") {
  
  
  message(paste("Running DL calling for",colnames(ctss)))
  ## Split on strand
  ctss.o <- ctss
  ctss <- calcPooled(ctss, inputAssay = callingAssay)
  ctss <- subset(ctss,score>0)
  
  #Call TCs for each sample
  object <- clusterUnidirectionally(ctss)
  object <- quantifyClusters(ctss,object)
  object <- calcPooled(object, inputAssay=inputAssay)
  object <- rowRanges(object) 
  
  message("Removing overlapping TCs by strand...")
  object <- swapRanges(object) #Summit focused
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
  
  
  # Summit mid calling 
  parallel_function <- function(c) {
    memb=con$membership[con$membership==c]
    d = object[names(memb)]
    db = d
    
    sm = d[strand(d)=="-"]
    sm = start(sm[sm$score==max(sm$score),])
    sm <- sm[1]  
    
    sp = d[strand(d)=="+"]
    sp = start(sp[sp$score==max(sp$score),])
    sp <- sp[length(sp)]  # If multiple summits, take the rightmost
    
    mid <- floor(sm + ((sp - sm) / 2))
    if (sm > sp) {  # if summits causes problem, revert to previous solution
      mid <- mergeLoci(end(d[strand(d) == "-"]), start(d[strand(d) == "+"]), d[strand(d) == "-"]$score, d[strand(d) == "+"]$score)
    }
    return(mid)
  }
  
  
  # Parallelize the function execution
  mid <- bplapply(unique(con$membership),parallel_function, BPPARAM = bpparam)
  
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

#k562.DLs.CTSS <- lapply(colnames(K562.r.CTSSs),function(x) divergentLociTCsSummit(K562.r.CTSSs[,x]))
k562.DLs.CTSS <- lapply(colnames(K562.r.CTSSs),function(x) divergentLociTCsSummit(K562.r.CTSSs[,x],callingAssay="counts.noSingletons"))

#K562
k562.s.CTSSs <- readRDS("../1.CTSSs/repSub.combined.k562.CTSSs.rds")

#Load in DHSs 
k562.dhs <- read.table("../0.External_resources/E123-DNase.macs2.hg38.narrowPeak")
colnames(k562.dhs) <- c("chr","start","end")
k562.dhs <- GRanges(k562.dhs[,1:3])
k562.dhs$thick <- IRanges(start=start(k562.dhs)+((end(k562.dhs)-start(k562.dhs))/2))
k562.dhs <- k562.dhs %>% swapRanges() %>% promoters(upstream=200, downstream=200) #Resize for equal width


#Annotating

olaps <- lapply(k562.DLs.CTSS, function(x) countOverlaps(x,k562.dhs)) #Count number of overlaps with DHSs for annotations
k562.DLs.CTSS <- suppressWarnings(lapply(k562.DLs.CTSS, function(x) assignTxType(x,txModels=txdb))) #Assign txType as annotations for later



#Some formating 
k562.DLs.CTSS <- GRangesList(k562.DLs.CTSS) 
k562.DLs.CTSS <- unlist(k562.DLs.CTSS)
k562.DLs.CTSS$dhs.olap <- unlist(olaps)


#Basic checks
tab <- table(k562.DLs.CTSS$sample, k562.DLs.CTSS$txType) %>% melt() %>% as.data.frame()
ggplot(tab, aes(Var1, value)) + geom_col() + coord_flip() + theme_bw() + scale_fill_manual(values=custom.col) + labs(x="sample", y="numberDLs")
ggplot(tab, aes(Var1, value,fill=Var2)) + geom_col(position="dodge") + coord_flip() + theme_bw() + scale_fill_manual(values=custom.col) + labs(x="sample", y="numberDLs")

#Tab stack
tab.stack <- tab %>%
  group_by(Var1) %>%
  reframe(total = sum(value),
            proportion = value / total,Var2=Var2) %>% 
  ungroup() 

ggplot(tab.stack, aes(Var1, proportion,fill=Var2)) + geom_col(position="stack") + coord_flip() + theme_bw() + scale_fill_manual(values=custom.col)+ labs(x="sample", y="proportions of DLs",fill="genomic annotations")




tab2 <- table(k562.DLs.CTSS$sample, k562.DLs.CTSS$dhs.olap) %>% melt() %>% as.data.frame()
ggplot(tab2, aes(Var1, value,fill=Var2>0)) + geom_col(position="dodge") + coord_flip() + theme_bw() + scale_fill_manual(values=custom.col) + labs(x="sample", y="numberDLs",fill="DHS overlap")


tab.stack2 <- tab2 %>%
  group_by(Var1) %>%
  reframe(total = sum(value),
          proportion = value / total,Var2=Var2) %>% 
  ungroup() 
ggplot(tab.stack2, aes(Var1, proportion,fill=Var2>0)) + geom_col(position="stack") + coord_flip() + theme_bw() + scale_fill_manual(values=custom.col)+ labs(x="sample", y="proportions of DLs",fill="DHS overlap")



k562.DLs.CTSS %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(directionality)) + geom_density() +facet_wrap(~sample) + theme_bw()


#Agreement between replicates?
seqlevels(k562.DLs.CTSS) <- seqlevels(K562.r.CTSSs)
k562.DLs <- quantifyDivergentLoci(k562.DLs.CTSS, K562.r.CTSSs, requireDisjoint = FALSE)

k562.DLs <- calcSupport(k562.DLs)

k562.DLs %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=txType)) + geom_histogram() + facet_wrap(~sample, scales="free_y") +scale_y_log10()



#Only nuclei
k562.DLs.n <- subset(k562.DLs[,colnames(k562.DLs) %in% c("K562_N4","K562_N5","K562_N6")], sample %in% c("K562_N4","K562_N5","K562_N6"))
k562.DLs.n <- calcSupport(k562.DLs.n)

k562.DLs.n %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=txType)) + geom_bar(position="dodge") + facet_wrap(~sample, scales="free_y") + theme_bw()+ scale_fill_manual(values=custom.col)


k562.DLs.n %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=dhs.olap>0)) + geom_bar(position="dodge") + facet_wrap(~sample, scales="free_y") + theme_bw()+ scale_fill_manual(values=custom.col)



#Only cellular
k562.DLs.c <- subset(k562.DLs[,colnames(k562.DLs) %in% c("K562_C1","K562_C2","K562_C3")], sample %in% c("K562_C1","K562_C2","K562_C3"))
k562.DLs.c <- calcSupport(k562.DLs.c)

k562.DLs.c %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=txType)) + geom_bar(position="dodge") + facet_wrap(~sample, scales="free_y") + theme_bw() + scale_fill_manual(values=custom.col)


k562.DLs.c %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=dhs.olap>0)) + geom_bar(position="dodge") + facet_wrap(~sample, scales="free_y") + theme_bw()+ scale_fill_manual(values=custom.col)

#Both set of replicates
k562.DLs.b <- subset(k562.DLs[,colnames(k562.DLs) %in% c("K562_C1","K562_C2","K562_C3","K562_N4","K562_N5","K562_N6")], sample %in% c("K562_C1","K562_C2","K562_C3","K562_N4","K562_N5","K562_N6","K562_groCap_1"))
k562.DLs.b <- calcSupport(k562.DLs.b)

k562.DLs.b %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=txType)) + geom_bar(position="dodge") + facet_wrap(~sample, scales="fixed") + theme_bw() + scale_fill_manual(values=custom.col)


k562.DLs.b %>% rowData() %>% as.data.frame(row.names=NULL) %>%
  ggplot(aes(support,fill=dhs.olap>0)) + geom_bar(position="dodge") + facet_wrap(~sample, scales="fixed") + theme_bw()+ scale_fill_manual(values=custom.col)





#Eqtls?
#Read in eQTL data from encode-rE2G 
#Downloaded from https://www.synapse.org/#!Synapse:syn52264240

bg.var <- import("../0.External_resources/all.bg.SNPs.hg38.baseline.v1.1.bed.sorted", format="bed")
bg.var <- assignTxType(bg.var, txModels = txdb)
bg.var <- subset(bg.var, txType%in% c("intron","intergenic"))

eqtls <- read.table(gzfile("../0.External_resources/GTEx_30tissues_release1.tsv.gz"),sep="\t")
colnames(eqtls) <- c("chromosome","start","end","hg19.ID","hg38.ID","ref","alt","cohort","method","tissue","gene","maf","beta_marginal","se_marginal","Zstats","pip","cs_id","beta_posterior","sd_posterior")

#Format to hg38
eqtls$end <- as.numeric(str_extract(eqtls$hg38.ID, "(?<=_)[0-9]+"))
eqtls$start <- eqtls$end-1 




eqtls05 <- subset(eqtls, pip>=0.5)
eqtls05 <- GRanges(eqtls05)
eqtls05 <- assignTxType(eqtls05,txModels = txdb)

eqtls05 <- subset(eqtls05,txType %in% c("intergenic","intron"))

eqtls05.grl <- split(eqtls05, eqtls05$tissue)

# Create a function to convert the eqtl data frame to a GRanges object
convert_to_GRanges <- function(df) {
  gr <- with(df, GRanges(seqnames = chromosome,
                         ranges = IRanges(start = start, end = end),
                         metadata = list(hg19.ID = hg19.ID,
                                         hg38.ID = hg38.ID,
                                         ref = ref,
                                         alt = alt,
                                         cohort = cohort,
                                         method = method,
                                         tissue = tissue,
                                         gene = gene,
                                         maf = maf,
                                         beta_marginal = beta_marginal,
                                         se_marginal = se_marginal,
                                         Zstats = Zstats,
                                         pip = pip,
                                         cs_id = cs_id,
                                         beta_posterior = beta_posterior,
                                         sd_posterior = sd_posterior)))
  return(gr)
}



#Split DL GRanges based on sample
k562.DLs.gr <- rowRanges(k562.DLs)
k562.DLs.gr <- keepStandardChromosomes(k562.DLs.gr,pruning.mode="coarse")
k562.DLs.gr <- subset(k562.DLs.gr,txType %in% c("intergenic"))
k562.DLs.grl <- split(k562.DLs.gr, mcols(k562.DLs.gr)$sample)
k562.DLs.grl <- lapply(k562.DLs.grl,function(gr){
  gr$rank = rank(-gr$score, ties.method = "random")
  gr = subset(gr, rank<1000)
  gr
})




# Function to count overlaps for each pair of GRanges
count_overlaps <- function(grl1, grl2) {
  overlaps <- lapply(seq_along(grl1), function(i) {
    lapply(seq_along(grl2), function(j) {
      # Find overlaps between GRanges in grl1 and grl2
      overlap <- findOverlaps(grl1[[i]], grl2[[j]], type = "any")
      return(length(overlap))
    })
  })
  return(overlaps)
}

# Count overlaps between the two GRangesList objects
overlaps <- count_overlaps(k562.DLs.grl , eqtls05.grl)


# Extract sample names from grl1
sample_names_grl1 <- names(k562.DLs.grl)

# Extract tissue names from grl2
tissue_names_grl2 <- names(eqtls05.grl)

# Create an empty data frame with row names from grl1 and column names from grl2
overlap_df <- matrix(unlist(overlaps), nrow = length(sample_names_grl1), 
                     byrow = TRUE, dimnames = list(sample_names_grl1, tissue_names_grl2))

# Convert to data frame
overlap_df <- as.data.frame(overlap_df)




#Calculate total number of eqtls
causal.total <- lapply(eqtls05.grl,length)
causal.total <- unlist(causal.total)

overlap_df_divided <- sapply(seq_along(causal.total), function(i) {
  overlap_df[, i] / causal.total[i]
  
})

overlap_df_divided = as.data.frame(overlap_df_divided)
row.names(overlap_df_divided) = row.names(overlap_df)
colnames(overlap_df_divided) = colnames(overlap_df)



#Add background calculations
bg.frac <- lapply(k562.DLs.grl , function(gr){
  overlap <- findOverlaps(gr, bg.var, type = "any")
  return(length(overlap)/length(bg.var))
})

overlap_df_divided_new <- sweep(overlap_df_divided, MARGIN = 1, STATS = unlist(bg.frac), FUN = "/")


k.full.pal <- c("#F5CE00", "#E2AA00", "#B68D00","#FF8D00","#DF613F","#DF5643","black","#92B9C9","#76A5B6","#004352","#00439B","#29A957","#008E37","#007A30","#4F5D5C")

overlap_df_divided_new$sample <- row.names(overlap_df_divided_new)
test <- melt(overlap_df_divided_new)
ggplot(test[test$sample%in%use,], aes(variable, value,color=sample)) + geom_point() + coord_flip() + theme_bw() + labs(y="enrichment", x="trait") + scale_color_manual(values=k.full.pal)


test2 <- data.frame(sample=overlap_df_divided_new$sample,enrichment=overlap_df_divided_new$Whole_Blood,recall=overlap_df_divided$Whole_Blood)

ggplot(test2[test2$sample %in%use,], aes(recall,enrichment,color=sample)) + geom_point() + scale_color_manual(values=k.full.pal) + theme_bw() + labs(x="recall",y="enrichment")
