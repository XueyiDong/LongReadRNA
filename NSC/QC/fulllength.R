source("Soneson/func.R")
library(ShortRead)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicFeatures)
library(Hmisc)

# read and calculate, save to RDS

dir_bam <- "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/results/rebasecall/minimap2_FLAMES_transcript"
bams <- paste0("barcode",c("07", 10:13, 15:17))

# make readDF
# for (i in 1:8){
#   cat("Reading", bams[i], "\n")
#   bam1 <- suppressWarnings(readBam(file.path(dir_bam, bams[i])))
#   cat("Making read DF\n")
#   readDF <- makeReadDf(bam1)
#   cat("Making summary list\n")
#   summ <- makeSummaryList(bam1)
#   saveRDS(readDF, file = file.path(dir_bam, paste0(bams[i], ".readDF.RDS")))
#   saveRDS(summ, file=file.path(dir_bam, paste0(bams[i], ".summ.RDS")))
# }

# read in RDS

readDF <- lapply(1:8, function(x){
  readRDS(file.path(dir_bam, paste0(bams[x], ".readDF.RDS")))})
DF <- rbind(readDF[[1]], readDF[[2]], readDF[[3]], readDF[[4]],
            readDF[[5]], readDF[[6]], readDF[[7]], readDF[[8]])
DF$sample <- rep(c("1", "2", "3", "4", "1", "5", "6", "7"),
                 c(nrow(readDF[[1]]), nrow(readDF[[2]]), nrow(readDF[[3]]),
                   nrow(readDF[[4]]), nrow(readDF[[5]]), nrow(readDF[[6]]),
                   nrow(readDF[[7]]), nrow(readDF[[8]])))
saveRDS(DF, file=file.path(dir_bam, "readDF.all.RDS"))
rm(readDF)




# read GFF and calculate tx len

# gff1 <- "/stornext/Home/data/allstaff/d/dong.x/analysis/2020/smchd1/NSC/flames/results/isoform_annotated.filtered.gff3"
# gff2 <- "/wehisan/home/allstaff/d/dong.x/annotation/Mouse/gencode.vM23.annotation.gff3"
# anno_FLAMES <- makeTxDbFromGFF(gff1, organism = "Mus musculus")
# anno_GENCODE <- makeTxDbFromGFF(gff2, organism = "Mus musculus")
# tl_FLAMES <- transcriptLengths(anno_FLAMES)
# tl_GENCODE <- transcriptLengths(anno_GENCODE)
# tl <- rbind(tl_FLAMES, tl_GENCODE)
# tl <- tl[!duplicated(tl$tx_name),]
# saveRDS(tl, "transcriptLengths.RDS")


# DF <- readRDS(file.path(dir_bam, "readDF.all.RDS"))
DF < -readRDS("../flames/results/readDF.RDS")
tl <- readRDS("./transcriptLengths.RDS")
m <- match(DF$seqnames, tl$tx_name)
DF$tx_len <- tl$tx_len[m]
# readDF <- na.omit(readDF)
# summary(readDF$width / readDF$tx_len)


# dataset info

# barcode <- read.delim("../scBarcode.tsv", stringsAsFactors = F)
# for(i in 1:nrow(barcode)){
#   tmp <- grepl(paste0("^", barcode$barcode[i]), readDF$read)
#   readDF$sample[tmp] <- barcode$name[i]
#   readDF$group[tmp] <- barcode$group[i]
# }

# make plots

maxLength = max(DF$tx_len)
DF$txLengthGroup <- cut2(DF$tx_len, cuts = c(0, 500, 1000, 1500, 
                                                     2000, maxLength))
DF$covFraction <- DF$width / DF$tx_len
DF$known <- grepl("^ENSMUST", DF$seqnames)
DF$group[DF$sample %in% c("1", "4", "5", "6")] <- "wt"
DF$group[DF$sample %in% c("2", "3", "7")] <- "Smchd1"


pdf("plots/fullLengthViolinLen.pdf")
# tx length
ggplot(DF, aes(x=txLengthGroup, y=covFraction)) +
  geom_violin() +
  theme_bw() +
  ylim(0, 1) 
dev.off()
pdf("plots/fullLengthViolinNovelty.pdf")
# tx length and novelty
ggplot(DF, aes(x=txLengthGroup, y=covFraction, fill = txLengthGroup)) +
  geom_violin() +
  theme_bw() +
  ylim(0, 1) 
dev.off()
pdf("fullLengthViolinSample.pdf")
# sample, filled by group
ggplot(DF, aes(x=sample, y=covFraction, fill=group)) +
  geom_violin() +
  theme_bw() +
  ylim(0, 1) +
  scale_fill_brewer(palette="Set2")
dev.off()
pdf("fullLengthBoxSample.pdf")
# sample, filled by group
ggplot(DF, aes(x=sample, y=covFraction, fill=group)) +
  geom_boxplot() +
  theme_bw() +
  ylim(0, 1) +
  scale_fill_brewer(palette="Set2")
dev.off()

