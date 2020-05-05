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


pdf("fullLengthViolinLen.pdf")
# tx length
ggplot(readDF, aes(x=txLengthGroup, y=covFraction, fill=txLengthGroup)) +
  geom_violin() +
  theme_bw() +
  ylim(0, 1) 
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
