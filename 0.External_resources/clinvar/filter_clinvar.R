library(readr)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cds <- reduce(unlist(cdsBy(txdb, "gene")))

files <- c(
  VUS="ClinVar_GRCh38_VUS_SNV.txt",
  LikelyPathogenic="ClinVar_GRCh38_LikelyPathogenic_SNV.txt",
  Pathogenic="ClinVar_GRCh38_Pathogenic_SNV.txt"
)

collapse_variants <- function(gr) {
  if (length(gr) == 0) return(gr)
  
  message("Collapsing ", length(gr), " variants into unique positions...")
  
  gr_unique <- gr
  mcols(gr_unique) <- NULL
  gr_unique <- unique(gr_unique)
  
  hits <- findOverlaps(gr_unique, gr)
  sh <- subjectHits(hits)
  qh <- queryHits(hits)
  
  ids_raw <- as.character(mcols(gr)$AlleleID)[sh]
  ids_split <- split(ids_raw, qh)
  mcols(gr_unique)$AlleleID <- unstrsplit(unique(CharacterList(ids_split)), sep = ";")
  
  sigs_raw <- as.character(mcols(gr)$ClinicalSignificance)[sh]
  sigs_split <- split(sigs_raw, qh)
  mcols(gr_unique)$ClinicalSignificance <- unstrsplit(unique(CharacterList(sigs_split)), sep = ";")
  
  if ("Category" %in% colnames(mcols(gr))) {
    cats_raw <- as.character(mcols(gr)$Category)[sh]
    cats_split <- split(cats_raw, qh)
    mcols(gr_unique)$Category <- unstrsplit(unique(CharacterList(cats_split)), sep = ";")
  }
  
  return(sort(gr_unique))
}

gr_list <- list()

for (cat in names(files)) {
  message("Processing ", cat)
  
  var_data <- read_tsv(files[[cat]], 
                       show_col_types=FALSE,
                       col_select = c("Chromosome", "Start", "Stop", "#AlleleID", "ClinicalSignificance"))
  
  var_data$Chromosome <- paste0("chr", var_data$Chromosome)
  
  gr <- makeGRangesFromDataFrame(
    var_data,
    seqnames.field="Chromosome",
    start.field="Start",
    end.field="Stop",
    keep.extra.columns=FALSE
  )
  
  start(gr) <- start(gr) + 1
  end(gr) <- end(gr) + 1
  
  gr$AlleleID <- as.character(var_data$`#AlleleID`)
  gr$ClinicalSignificance <- as.character(var_data$ClinicalSignificance)
  gr$Category <- as.character(cat) # Need this for the final merged collapse
  
  rm(var_data)
  gc()
  
  gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
  
  olap <- findOverlaps(gr, cds)
  if (length(olap) > 0) {
    gr <- gr[-queryHits(olap)]
  }
  
  gr <- collapse_variants(gr)
  
  gr_out <- gr
  mcols(gr_out)$Category <- NULL
  
  prefix <- paste0("ClinVar_GRCh38_", cat, "_noncoding_SNV")
  write_tsv(as.data.frame(gr_out), paste0(prefix, ".tsv"))
  export(gr_out, paste0(prefix, ".bed"), format="BED")
  
  gr_list[[cat]] <- gr
}

gr_list <- Filter(function(x) length(x) > 0, gr_list)
gr_all <- do.call(c, unname(gr_list))
gr_all <- collapse_variants(gr_all)

gr_all_out <- gr_all
mcols(gr_all_out)$Category <- NULL

write_tsv(as.data.frame(gr_all_out), "ClinVar_GRCh38_allCategories_noncoding_SNV.tsv")
export(gr_all_out, "ClinVar_GRCh38_allCategories_noncoding_SNV.bed", format="BED")