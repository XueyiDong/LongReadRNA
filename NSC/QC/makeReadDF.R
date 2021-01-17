source("func.R")
library(ShortRead)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicFeatures)
library(Hmisc)

dir <- "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/results/rebasecall/minimap2_FLAMES_transcript"
samples <- paste0("barcode", c("07", 10:13, 15:17))
bams <- paste0(samples, ".sorted.bam")

for (i in 1:8){
  cat("Reading bam file.\n")
  bam1 <- suppressWarnings(readBam(file.path(dir, bams[i])))
  cat("Making read DF.\n")
  readDF <- makeReadDf(bam1)
  cat("Making summary list.\n")
  summ <- makeSummaryList(bam1)
  cat("Saving readDF.\n")
  saveRDS(readDF, paste0(dir, "/", samples[i], ".readDF.RDS"))
  cat("Saving summary list.\n")
  saveRDS(summ, paste0(dir, "/", samples[i], ".summList.RDS"))
}