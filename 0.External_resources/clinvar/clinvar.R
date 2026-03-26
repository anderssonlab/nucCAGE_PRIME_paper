
library(readr)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(tidyverse)

var_data <- read_tsv("ClinVar_GRCh38_VUS_SNV.txt")
var_data$Chromosome <- paste0("chr",var_data$Chromosome)

gr <- makeGRangesFromDataFrame(var_data)
gr$AlleleID <- var_data$`#AlleleID`
gr$Name <- var_data$Name
gr$ClinicalSignificance <- var_data$ClinicalSignificance
start(gr) <- start(gr) + 1
end(gr) <- end(gr) + 1
gr <- keepStandardChromosomes(gr, pruning.mode="coarse")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cds <- unlist(cdsBy(txdb,"gene"))
olap <- findOverlaps(gr, cds, minoverlap=1)
gr <- gr[-queryHits(olap)]
gr <- sort(gr)

olap <- countOverlaps(gr,gr)
new.gr <- gr[which(olap==1)]

nonmerged.gr <- gr[which(olap>1)]
pos <- paste0(seqnames(nonmerged.gr),":",start(nonmerged.gr))
pos.unique <- unique(pos)
df <- as.data.frame(nonmerged.gr)[,c("AlleleID","Name","ClinicalSignificance")]
df <- do.call("rbind",by(df, pos, function(x) sapply(x,paste,collapse=";")))

merged.gr <- nonmerged.gr[match(pos.unique,pos)]
merged.gr$AlleleID <- df[,"AlleleID"]
merged.gr$Name <- df[,"Name"]
merged.gr$ClinicalSignificance <- df[,"ClinicalSignificance"]

new.gr <- sort(c(new.gr,merged.gr))

write_tsv(as.data.frame(new.gr),"ClinVar_GRCh38_VUS_SNV.tsv")

ClinVar.data <- new.gr

col.scheme <- base::c(
  "#377eb8","#d39200","#93aa00",
           "#4daf4a","#984ea3","#ff7f01",
           "#e41b1c","#256866","#f8bc24",
           "grey","black","hotpink"
)

## Read in SNP.file
SNP.data <- rtracklayer::import.bed( # 9.991.229
  con = base::paste0(
    here::here("/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/"),
    "all.bg.SNPs.hg38.baseline.v1.1.bed.sorted"
  )
)
SNP.data <- BiocGenerics::sort(SNP.data)
olap <- GenomicRanges::findOverlaps( ## 93.742
  query = SNP.data,
  subject = cds,
  minoverlap = 1L
)
SNP.data$name <- NULL
SNP.data <- base::unique( # 9.914.129
  SNP.data[-S4Vectors::queryHits(olap)]
)
GenomeInfoDb::seqlevels(txdb, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(SNP.data)

#ENCODE cCREs
cCRE <- readr::read_tsv(
  file = base::paste0("",
    "GRCh38-cCREs.bed"
  ),
  col_names = FALSE,
  show_col_types = FALSE
)
df <- cCRE[,1:3]
base::colnames(df) <- base::c("chrom","start","end")
cCRE.data <- GenomicRanges::makeGRangesFromDataFrame(df)
cCRE.data$class <- cCRE[,6, drop = TRUE]

promoters <- IRanges::promoters(
  x = txdb,
  upstream = 500,
  downstream = 500
)

point.data <- list()

for (n in base::c("dELS","pELS","PLS")) {
  x <- base::subset(cCRE.data, class == n)
  x$class <- NULL
  GenomeInfoDb::seqlevels(x,pruning.mode="coarse") <- GenomeInfoDb::seqlevels(SNP.data)
  point.data[[paste0("ENCODE cCRE ",n)]] <- x
}

point.data[["ENCODE cCRE"]] <- base::c(
  point.data[["ENCODE cCRE dELS"]],
  point.data[["ENCODE cCRE pELS"]],
  point.data[["ENCODE cCRE PLS"]]
)

PRIME.proximal <- readr::read_tsv(
  file = base::paste0("/maps/projects/ralab/data/projects/nucleiCAGEproject/resource/PRIME_FANTOM5_agnostic/",
                      "PRIME_FANTOM5_agnostic_proximal_0.5.bed"
  ),
  col_names = FALSE,
  show_col_types = FALSE
)
PRIME.proximal <- PRIME.proximal[,1:3]
base::colnames(PRIME.proximal) <- base::c("chrom","start","end")
PRIME.proximal <- GenomicRanges::makeGRangesFromDataFrame(PRIME.proximal)
point.data[["PRIME proximal"]] <- PRIME.proximal

PRIME.distal <- readr::read_tsv(
  file = base::paste0("/maps/projects/ralab/data/projects/nucleiCAGEproject/resource/PRIME_FANTOM5_agnostic/",
                      "PRIME_FANTOM5_agnostic_distal_0.5.bed"
  ),
  col_names = FALSE,
  show_col_types = FALSE
)
PRIME.distal <- PRIME.distal[,1:3]
base::colnames(PRIME.distal) <- base::c("chrom","start","end")
PRIME.distal <- GenomicRanges::makeGRangesFromDataFrame(PRIME.distal)
point.data[["PRIME distal"]] <- PRIME.distal

point.data[["PRIME"]] <- c(point.data[["PRIME distal"]],
                           point.data[["PRIME proximal"]])

resize.width = 200
for (n in base::names(point.data)) {
  point.data[[n]] <- base::sort(point.data[[n]])
  GenomeInfoDb::seqlevels(point.data[[n]],pruning.mode="coarse") <- GenomeInfoDb::seqlevels(SNP.data)
  point.data[[n]] <- IRanges::resize(
    x = point.data[[n]],
    width = resize.width,
    fix = "center"
  )
  point.data[[n]] <- IRanges::reduce(
    x = point.data[[n]],
    min.gapwidth = 0
  )
}

point.enrichment.precision.recall <- function(data, fg.data, bg.data) {
  
  base::do.call(
    what = "rbind",
    args = base::lapply(
      X = base::names(data),
      FUN = function(n) {
        
        fg.olap <- GenomicRanges::findOverlaps(
          query = fg.data,
          subject = data[[n]],
          minoverlap = 1
        )
        fg.count <- base::length(
          base::unique(
            S4Vectors::queryHits(fg.olap)
          )
        )
        region.count <- base::length(
          base::unique(
            S4Vectors::subjectHits(fg.olap)
          )
        )
        bg.olap <- GenomicRanges::findOverlaps(
          query = bg.data,
          subject = data[[n]],
          minoverlap = 1
        )
        bg.count <- base::length(
          base::unique(
            S4Vectors::queryHits(bg.olap)
          )
        )
        return(
          base::data.frame(
            enrichment = (fg.count/base::length(fg.data)) / (bg.count/base::length(bg.data)),
            precision = fg.count/(fg.count + bg.count),
            recall = fg.count/base::length(fg.data),
            fg.count = fg.count,
            bg.count = bg.count,
            region.count = region.count,
            sample = n
          )
        )
      }
    )
  )
}

point.ClinVar.df <- point.enrichment.precision.recall(
  data = point.data,
  fg.data = ClinVar.data,
  bg.data = SNP.data
)

VUS.n <- length(ClinVar.data)

ggplot2::ggplot(
  data = point.ClinVar.df,
  mapping = aes(
    x = recall,
    y = enrichment,
    color = sample,
    label = sample
  )
) + 
  geom_point() +
  theme_bw() +
  theme(
    aspect.ratio = 1
  ) +
  ggrepel::geom_text_repel(
    size = 2,
    min.segment.length = 0,
    color = "black"
  ) + 
  xlab(
    base::paste0("recall (n=",VUS.n,")")
  ) +
  scale_color_manual(
    values = col.scheme
  )

ggplot2::ggsave(
  filename = base::paste0(
    here::here(""),
    "VUS_enrichment_vs_recall_point_estimates.pdf"
  ),
  width = 5,
  height = 3,
  device = "pdf"
)

ggplot2::ggplot(
  data = subset(point.ClinVar.df,sample %in% c("ENCODE cCRE dELS","PRIME distal")),
  mapping = aes(
    x = recall,
    y = enrichment,
    color = sample,
    label = sample
  )
) + 
  geom_point() +
  theme_bw() +
  theme(
    aspect.ratio = 1
  ) +
  ggrepel::geom_text_repel(
    size = 2,
    min.segment.length = 0,
    color = "black"
  ) + 
  xlab(
    base::paste0("recall (n=",VUS.n,")")
  ) +
  scale_color_manual(
    values = col.scheme
  )

ggplot2::ggsave(
  filename = base::paste0(
    here::here(""),
    "VUS_enrichment_vs_recall_point_estimates_distal.pdf"
  ),
  width = 5,
  height = 3,
  device = "pdf"
)
