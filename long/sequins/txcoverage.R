source("../run2/QC/Soneson/func.R")
library(ShortRead)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicFeatures)
library(Hmisc)

# read and calculate, save to RDS

dir_bam <- "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/realign2transcript.bam"

bam1 <- suppressWarnings(readBam(dir_bam))
readDF <- makeReadDf(bam1)
summ <- makeSummaryList(bam1)
saveRDS(readDF, file = "sequins.readDF.RDS")
saveRDS(summ, file= "sequins.summ.RDS")
saveRDS(bam1, file= "bam.RDS")
rm(bam1)
readDF <- readRDS("sequins.readDF.RDS")
# attach tx len from annotation to readDF

readDF$tx_len <- anno$LENGTH[match(readDF$seqnames, anno$NAME)]

# remove 'novel' isoforms
readDF$novelty <- is.na(readDF$tx_len)
readDF <- readDF[!readDF$novelty,]


maxLength = max(readDF$tx_len)
readDF$txLengthGroup <- cut2(readDF$tx_len, cuts = c(0, 500, 1000, 1500, 
                                             2000, maxLength))
readDF$covFraction <- readDF$width / readDF$tx_len
# DF$group[DF$sample %in% c("1", "4", "5", "6")] <- "wt"
# DF$group[DF$sample %in% c("2", "3", "7")] <- "Smchd1"


pdf("fullLengthViolinLen.pdf", height = 5, width = 8)
give.n <- function(x){
  return(c(y = 0.98, label = length(x))) 
}
# tx length
ggplot(readDF, aes(x=txLengthGroup, y=covFraction, fill=txLengthGroup)) +
  geom_violin() +
  theme_bw() +
  ylim(0, 1) +
  stat_summary(fun.data = give.n, geom = "text", vjust = -1)
dev.off()

pdf("fullLengthBoxLen.pdf", height = 5, width = 8)
give.n <- function(x){
  return(c(y = 0.98, label = length(x))) 
}
# tx length
ggplot(readDF, aes(x=txLengthGroup, y=covFraction, fill=txLengthGroup)) +
  geom_boxplot() +
  theme_bw() +
  ylim(0, 1) +
  stat_summary(fun.data = give.n, geom = "text", vjust = -1)
dev.off()

library(viridis)
pdf("mappedVsTxLen.pdf")
smoothScatter(readDF$tx_len, readDF$width, 
              xlab = "transcript length",
              ylab = "mapped length",
              log = "xy")
abline(coef = c(0,1), col="red")
# ggplot(readDF, aes(x=tx_len, y=width)) +
#   geom_hex(binwidth = c(0.1, 0.1)) +
#   theme_bw() +
#   labs(x = "transcript length", y = "mapped length") +
#   scale_fill_viridis_c(direction = -1) +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10') +
#   geom_abline(intercept = 0, slope = 1, colour="red")
dev.off()

# thought: normal scatter plot, colour by cov frection. Problem: too many dots. may need some single cell package to solve.

# look into why some cov fraction > 1

covFracLargerThan1 <- DF[DF$covFraction>1, ]
head(covFracLargerThan1)
dim(covFracLargerThan1)
length(unique(covFracLargerThan1$read))
length(unique(covFracLargerThan1$seqnames))

# make alternative plot
pdf("txLenCovFrac.pdf", height = 5, width = 8)
ggplot(readDF, aes(x=tx_len, y=covFraction)) +
  stat_binhex(binwidth = c(0.03, 0.1)) +
  theme_bw() +
  labs(x="Transcript length", y="Coverage fraction") +
  scale_x_continuous(trans = "log10") +
  scale_fill_viridis(direction = -1, option="A")
dev.off()

#calculate some stats by transcript

txStat <- sapply(unique(readDF$seqnames), function(x){
  readDF.sel = readDF[readDF$seqnames==x, ]
  meanCovFrac = mean(readDF.sel$covFraction)
  medianCovFraction = median(readDF.sel$covFraction)
  fl95 = sum(readDF.sel$covFraction >= 0.95) / nrow(readDF.sel)
  fl90 = sum(readDF.sel$covFraction >= 0.90) / nrow(readDF.sel)
  c(readDF.sel[1, "tx_len"], meanCovFrac, medianCovFraction, fl95, fl90, nrow(readDF.sel))
}, simplify = TRUE)

txStat <- as.data.frame(t(txStat))
colnames(txStat) <- c("tx_len", "mean", "median", "fl95", "fl90", "count")
txStat$log_count <- log(txStat$count)

ggplot(txStat, aes(x=tx_len, y=median, size=log_count))+
  scale_x_continuous(trans = "log10") +
  geom_point()

ggplot(txStat, aes(x=tx_len, y=mean, size=log_count))+
  scale_x_continuous(trans = "log10") +
  geom_point()

ggplot(txStat, aes(x=tx_len, y=fl90, size=log_count))+
  scale_x_continuous(trans = "log10") +
  geom_point()

pdf("txLenFL.pdf", height = 5, width = 8)
ggplot(txStat, aes(x=tx_len, y=fl95, colour = log_count))+
  scale_x_continuous(trans = "log10") +
  geom_point() +
  labs(x = "Annotated transcript length", y = "Fraction of full-length", colour = "log count") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_colour_viridis()
dev.off()


readDF$alignedFraction <- readDF$alignedLength / readDF$readLength

pdf("alignCov.pdf")
ggplot(readDF, aes(x=alignedFraction, y=covFraction)) +
  stat_binhex()+
  scale_fill_viridis(direction = -1, option="A")
dev.off()

smoothScatter(readDF$alignedFraction, readDF$covFraction)
cor(readDF$alignedFraction, readDF$tx_len)
#-0.01392428 almost no corr between tx len and aligned fraction

cor(readDF$tx_len, readDF$covFraction)
# -0.5247596

cor(readDF$alignedFraction, readDF$covFraction)
# 0.6427299


smoothScatter(readDF$tx_len, readDF$covFraction, log="x")

readDF$isFullLength <- readDF$covFraction >= 0.95

pdf("AlignFL.pdf", height = 5, width = 8)
ggplot(readDF, aes(x=isFullLength, y=alignedFraction, fill=isFullLength)) +
  geom_violin() +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "none") +
  labs(x = "Is full-length read", y = "Fraction of aligned bases")
dev.off()

