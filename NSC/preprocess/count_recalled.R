# Sample information
targets <- as.data.frame(cbind(barcode=7:18, group=c("wt", "fresia", "fresia", "smchd1", "smchd1", "wt", "wt", "fresia", "wt", "wt", "smchd1", "fresia"), batch=rep(1:2, each=6)))
# targets$sample <- paste0(targets$group, c(1,1,2,1,2,2))

# Set main
# setwd("/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2/")

# Bam files
path1.bam <- "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/Xueyi/results/minimap2_mm10/"
path2.bam <- "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2/results/minimap2_mm10/"
files.bam <- c(paste0("BC0", targets$barcode[1:3], ".sorted.bam"), paste0("BC", targets$barcode[4:12], ".sorted.bam"))
files <- c(paste0(path1.bam, files.bam[1:6]), paste0(path2.bam, files.bam[7:12]))

#-------------------------------------- fc
library(Rsubread)
fc <- featureCounts(files=files, annot.inbuilt="mm10", allowMultiOverlap=FALSE, nthreads=8, isLongRead=TRUE,primaryOnly=TRUE)
save(fc, file="countsPrimary.RData")
colnames(fc$stat)=files.bam
sum(fc$stat[1,-1])
# colnames(c2$stat) <- c("Status", paste0("BC", targets$barcode))
# c2$stat

# save(c1,c2, file = "~/analysis/smchd1/long/run2/counts.RData")

#--------------------------------------- use GenomicAlignments
# read in bam files
library(GenomicAlignments)
filepath <- paste0(path.bam, files.bam)
aln <- lapply(filepath, readGAlignments)
# read in features
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
GR <- exonsBy(txdb, by="gene")
# summarize overlaps
olap <- lapply(aln, function(x){
  summarizeOverlaps(GR, x, ignore.strand = T)
})
# organize results
counts <- data.frame(sapply(olap, assay, simplify = T))
rownames(counts) <- rownames(assay(olap[[1]]))
colnames(counts) <- paste0("BC", targets$barcode)
# save results
setwd("~/analysis/smchd1/long/run2")
save(olap, file="summarizeOverlaps.RData")
save(counts, file="counts_olap.RData")
