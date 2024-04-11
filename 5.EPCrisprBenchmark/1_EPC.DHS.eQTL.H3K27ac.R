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
  
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

## 0. Setting
Dis.EPC <- 250
Dis.DHS <- 200
dir.resources <- "/projects/ralab/data/projects/nucleiCAGEproject/0.K562_EPCrisprBenchmark_resources"

genome.info <- keepStandardChromosomes(SeqinfoForUCSCGenome("hg38"), species="Homo sapiens")
standard.chromosomes <- seqnames(genome.info)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#### EPC 

# 1. import file, filter out NA, mutate (start, end), and keep std chr
EPC.raw <- read_tsv(paste0(dir.resources, "/EPCrisprBenchmark_ensemble_data_GRCh38.tsv"), show_col_types = FALSE)
EPC.df <- EPC.raw %>% 
  dplyr::select(c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "Regulated", "EffectSize", "Significant", "pValueAdjusted")) %>% 
  dplyr::filter(Regulated != "NA", Significant != "NA") 

colnames(EPC.df) <- c("seqnames", "start", "end", "measuredGeneSymbol", "Regulated", "EffectSize", "Significant", "pValueAdjusted")

EPC.df <- EPC.df %>%
  dplyr::mutate(thick = start + round((end - start + 1) / 2),
                end = thick + Dis.EPC,
                start = thick - Dis.EPC,
                start = dplyr::case_when(start < 0 ~ 0,
                                         TRUE ~ start)) %>%
  dplyr::filter(seqnames %in% standard.chromosomes)

# 2. categorize data based on regulate and significant
pos.df <- EPC.df %>% dplyr::filter(Regulated == TRUE)
neg.notsig.df <- EPC.df %>% dplyr::filter(Regulated == FALSE, Significant == FALSE)
neg.sigfcpos.df <- EPC.df %>% dplyr::filter(Regulated == FALSE, Significant == TRUE, EffectSize >= 0)

# 3. grouping the 'exact' same identified enhancers 
pos.gr <- pos.df %>% group_by(seqnames, thick) %>% summarise_all(funs(paste(unique(.), collapse = ";"))) %>% GRanges()
neg.notsig.gr <- neg.notsig.df %>% group_by(seqnames, thick) %>% summarise_all(funs(paste(unique(.), collapse = ";"))) %>% GRanges()
neg.sigfcpos.gr <- neg.sigfcpos.df %>% group_by(seqnames, thick) %>% summarise_all(funs(paste(unique(.), collapse = ";"))) %>% GRanges()

# 4. remove overlap 
# # prioritization: pos > neg_sigfcpos > neg_notsig
neg.notsig.gr <- subsetByOverlaps(neg.notsig.gr, c(pos.gr, neg.sigfcpos.gr), type="any", invert=TRUE) 

# 5. annotation
# previously set to tssUpstream=500, tssDownstream=500
# default was set (tssUpstream = 100, tssDownstream = 100) to match hjolli's
pos.gr <- CAGEfightR::assignTxType(pos.gr, txModels=txdb) 
neg.notsig.gr <- CAGEfightR::assignTxType(neg.notsig.gr, txModels=txdb)

pos.gr.txType <- table(pos.gr$txType) / length(pos.gr) * 100 
pos.plt <- data.frame(pos.gr.txType)
colnames(pos.plt) <- c("txType", "percentage")

neg.notsig.gr.txType <- table(neg.notsig.gr$txType) / length(neg.notsig.gr) * 100
neg.plt <- data.frame(neg.notsig.gr.txType)
colnames(neg.plt) <- c("txType", "percentage")

combined.plt <- rbind(
  transform(pos.plt, group = paste0("Positive (total count = ", as.character(length(pos.gr)), ")")),
  transform(neg.plt, group = paste0("Negative (total count = ", as.character(length(neg.notsig.gr)), ")")))

p <- ggplot(combined.plt, aes(x = txType, y = percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage of EPCrisprBenchmark txType",
       x = "txType",
       y = "percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("1.Percentage_of_EPCrisprBenchmark_txType.pdf", plot = p, width = 8, height = 6, units = "in")

# 6. save .RData for profiles
saveRDS(pos.gr, "1.EPC.pos.gr.rds")
saveRDS(neg.notsig.gr, "1.EPC.neg.notsig.gr.rds")



#### DHS

# 1. prep DHS
DHS.raw <- read_tsv(paste0(dir.resources, "/E123-DNase.macs2.hg38.narrowPeak"), 
                   col_names = c("seqnames","start","end","rank"),
                   col_types = readr::cols(.default = col_double(), seqnames = col_character(), rank = col_character())) 
DHS.df <- DHS.raw %>% data.frame() %>% 
  dplyr::mutate(thick = start + round((end - start + 1) / 2),
                end_new = thick + Dis.DHS,
                start_new = thick - Dis.DHS,
                start_new = dplyr::case_when(start_new < 0 ~ 0,
                                           TRUE ~ start_new)) %>%
  dplyr::select(seqnames, start = "start_new", end = "end_new", thick) %>%
  dplyr::filter(seqnames %in% standard.chromosomes)
DHS.gr <- GRanges(DHS.df)

# 2. overlap DHS
pos.gr$count.DHS = countOverlaps(pos.gr, DHS.gr)
neg.notsig.gr$count.DHS = countOverlaps(neg.notsig.gr, DHS.gr)

# 3. decision making
pos.gr$ovl.DHS = ifelse(pos.gr$count.DHS > 0, TRUE, FALSE)
neg.notsig.gr$ovl.DHS = ifelse(neg.notsig.gr$count.DHS > 0, TRUE, FALSE)



#### fine-mapped eQTL

# 1. prep eQTL
eQTL.raw <- read_tsv(paste0(dir.resources, "/eQTL/QTD000549.credible_sets.tsv"), show_col_types = FALSE)
eQT.df <- eQTL.raw  %>% data.frame() %>% 
  dplyr::select(variant, pip) %>% 
  separate(col = variant, into = c("seqnames", "start", "from", "to"), sep = "_") %>%
  dplyr::mutate(end=start)  %>% 
  dplyr::select(seqnames, start, end, pip)

eQTL.neg0.01.gr <- eQT.df %>% dplyr::filter(pip < 0.01) %>% GRanges()
eQTL.pos0.5.gr <- eQT.df %>% dplyr::filter(pip > 0.5) %>% GRanges()
eQTL.pos0.9.gr <- eQT.df %>% dplyr::filter(pip > 0.9) %>% GRanges()

# 2. count overlap eQTL
pos.gr$count.eQTL.neg0.01 = countOverlaps(pos.gr, eQTL.neg0.01.gr)
neg.notsig.gr$count.eQTL.neg0.01 = countOverlaps(neg.notsig.gr, eQTL.neg0.01.gr)

pos.gr$count.eQTL.pos0.5 = countOverlaps(pos.gr, eQTL.pos0.5.gr)
neg.notsig.gr$count.eQTL.pos0.5 = countOverlaps(neg.notsig.gr, eQTL.pos0.5.gr)

pos.gr$count.eQTL.pos0.9 = countOverlaps(pos.gr, eQTL.pos0.9.gr)
neg.notsig.gr$count.eQTL.pos0.9 = countOverlaps(neg.notsig.gr, eQTL.pos0.9.gr)

# 3. decision making

pos.gr$ovl.eQTL0.5 <- ifelse(pos.gr$count.eQTL.pos0.5 > 0, TRUE,
                             ifelse(pos.gr$count.eQTL.neg0.01 > 0, FALSE, NA))
neg.notsig.gr$ovl.eQTL0.5 <- ifelse(neg.notsig.gr$count.eQTL.pos0.5 > 0, TRUE,
                                    ifelse(neg.notsig.gr$count.eQTL.neg0.01 > 0, FALSE, NA))

pos.gr$ovl.eQTL0.9 <- ifelse(pos.gr$count.eQTL.pos0.9 > 0, TRUE,
                             ifelse(pos.gr$count.eQTL.neg0.01 > 0, FALSE, NA))
neg.notsig.gr$ovl.eQTL0.9 <- ifelse(neg.notsig.gr$count.eQTL.pos0.9 > 0, TRUE,
                                    ifelse(neg.notsig.gr$count.eQTL.neg0.01 > 0, FALSE, NA))


#### save last object
#save(list=c("pos.gr", "neg.notsig.gr"), file=paste0("1.EPC.DHS.eQTL.RData"))

#write.csv(data.frame(pos.gr), "1.EPC.DHS.eQTL.pos.csv")
#write.csv(data.frame(neg.notsig.gr), "1.EPC.DHS.eQTL.neg.csv")


#### H3K27ac
#library(rtracklayer)
#K27ac.raw <- rtracklayer::import.bw(paste0(dir.resources, "/H3K27ac_k562_hg38/J599.ChIP-seq_H3K27ac_signal_p-value.bigwig"))
#bw.gr=import(bwFile, which=promoter.gr)



