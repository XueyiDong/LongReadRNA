# Sample information
# targets <- as.data.frame(cbind(barcode=7:18, group=c("wt", "fresia", "fresia", "smchd1", "smchd1", "wt", "wt", "fresia", "wt", "wt", "smchd1", "fresia"), batch=rep(1:2, each=6)))
# targets$sample <- paste0(targets$group, c(1,1,2,1,2,2))

# Set main
# setwd("/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA_rep2/")

# Bam files
path <- "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/results/rebasecall/minimap2_mm10"
files <-file.path(path, list.files(path))
anno <- "/wehisan/home/allstaff/d/dong.x/annotation/Mouse/gencode.vM23.annotation.gff3"

#-------------------------------------- fc
library(Rsubread)
fc <- featureCounts(files=files, 
                    annot.ext = anno,
                    isGTFAnnotationFile = TRUE,
                    allowMultiOverlap=FALSE, 
                    nthreads=8, 
                    isLongRead=TRUE,
                    primaryOnly=TRUE)
saveRDS(fc, file="counts.RDS")

 sum(fc$stat[1,-1])
# colnames(c2$stat) <- c("Status", paste0("BC", targets$barcode))
# c2$stat

# save(c1,c2, file = "~/analysis/smchd1/long/run2/counts.RData")


