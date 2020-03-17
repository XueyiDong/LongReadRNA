source("Soneson/func.R")
library(ShortRead)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicFeatures)
library(Hmisc)

# read and calculate, save to RDS

dir_bam <- "/wehisan/general/old-prkfs2/disk503/GP_Transfer/Smchd1/long/scFLT/run1"
bams <- paste0("Tx.BC", c("07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18"), ".bam")
for (i in 1:12){
  cat("Reading", bams[i], "\n")
  bam1 <- suppressWarnings(readBam(file.path(dir_bam, bams[i])))
  cat("Making read DF\n")
  readDF <- makeReadDf(bam1)
  cat("Making summary list\n")
  summ <- makeSummaryList(bam1)
  saveRDS(readDF, file = file.path(dir_bam, paste0(bams[i], "readDF.RDS")))
  saveRDS(summ, file=file.path(dir_bam, paste0(bams[i], ".summ.RDS")))
}

# read in RDS

readDF <- lapply(c(1, 4, 5, 6, 7, 9, 10, 11), function(x){
  readRDS(file.path(dir_bam, paste0(bams[x], "readDF.RDS")))})
DF <- rbind(readDF[[1]], readDF[[2]], readDF[[3]], readDF[[4]],
            readDF[[5]], readDF[[6]], readDF[[7]], readDF[[8]])
DF$sample <- rep(c(1, 2, 3, 4, 1, 5, 6, 7), 
                 c(nrow(readDF[[1]]), nrow(readDF[[2]]), nrow(readDF[[3]]),
                   nrow(readDF[[4]]), nrow(readDF[[5]]), nrow(readDF[[6]]),
                   nrow(readDF[[7]]), nrow(readDF[[8]])))
saveRDS(DF, file=file.path(dir_bam, "readDF.all.RDS"))
rm(readDF)


# read GFF and calculate tx len

gff1 <- "/wehisan/general/old-prkfs2/disk503/GP_Transfer/Smchd1/long/scFLT/run1/isoform_annotated.gff3"
gff2 <- "/wehisan/home/allstaff/d/dong.x/annotation/Mouse/gencode.vM23.annotation.gff3"
anno_FLTSA <- makeTxDbFromGFF(gff1, organism = "Mus musculus")
anno_GENCODE <- makeTxDbFromGFF(gff2, organism = "Mus musculus")
tl_FLTSA <- transcriptLengths(anno_FLTSA)
tl_GENCODE <- transcriptLengths(anno_GENCODE)
tl <- rbind(tl_FLTSA, tl_GENCODE)
tl <- tl[!duplicated(tl$tx_name),]
save(tl, file = "/wehisan/general/old-prkfs2/disk503/GP_Transfer/Smchd1/long/scFLT/run1//rdata/tl.RData")

load("/wehisan/general/old-prkfs2/disk503/GP_Transfer/Smchd1/long/scFLT/run1/rdata/tl.RData")


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

maxLength = max(readDF$tx_len)
readDF$txLengthGroup <- cut2(readDF$tx_len, cuts = c(0, 500, 1000, 1500, 
                                                     2000, maxLength))
readDF$covFraction <- readDF$width / readDF$tx_len
readDF$known <- grepl("^ENSMUST", readDF$seqnames)



pdf("fullLengthViolin.pdf")
# tx length
ggplot(readDF, aes(x=txLengthGroup, y=covFraction)) +
  geom_violin() +
  theme_bw()
# tx length and novelty
ggplot(readDF, aes(x=txLengthGroup, y=covFraction, fill = txLengthGroup)) +
  geom_violin() +
  theme_bw()
# sample, fill by group
ggplot(readDF, aes(x=sample, y=covFraction, fill=group)) +
  geom_violin() +
  theme_bw()
# group, fill by group
ggplot(readDF, aes(x=group, y=covFraction, fill=group)) +
  geom_violin() +
  theme_bw()
dev.off()

# plots wo fresia

readDF <- readDF[readDF$group == "wt" | readDF$group == "smchd1",]
save(readDF, file = "readDF_filtered.RData")

pdf("fullLengthViolinNoFr.pdf")
# tx length
ggplot(readDF, aes(x=txLengthGroup, y=covFraction)) +
  geom_violin() +
  theme_bw()
# tx length and novelty
ggplot(readDF, aes(x=txLengthGroup, y=covFraction, fill = txLengthGroup)) +
  geom_violin() +
  theme_bw()
# sample, fill by group
ggplot(readDF, aes(x=sample, y=covFraction, fill=group)) +
  geom_violin() +
  theme_bw()
# group, fill by group
ggplot(readDF, aes(x=group, y=covFraction, fill=group)) +
  geom_violin() +
  theme_bw()
dev.off()


