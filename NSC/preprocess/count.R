# Bam files
path <- "/wehisan/general/academic/seq_data/quentin/Nanopore/Smchd1-NSC-cDNA/results/rebasecall/minimap2_mm10"
files <-file.path(path, list.files(path))
anno <- "/wehisan/home/allstaff/d/dong.x/annotation/Mouse/gencode.vM23.annotation.gff3"

library(Rsubread)
fc <- featureCounts(files=files, 
                    annot.ext = anno,
                    isGTFAnnotationFile = TRUE,
                    allowMultiOverlap=FALSE, 
                    nthreads=8, 
                    isLongRead=TRUE,
                    primaryOnly=TRUE)
saveRDS(fc, file="counts.RDS")