dir <- "/home/users/allstaff/mritchie/mattlab/long_RNA_benchmark/Illumina_April19"
libs <- c(paste(1:2, "-H1975-", 1:2, "_S", 1:2, sep=""), paste(6:7, "-HCC-", 1:2, "_S", 6:7, sep=""))
bam <- paste(dir, "/results/subread_sequins/", libs, ".sorted.bam", sep="")
gtf <- "~/annotation/sequins/rnasequin_annotation_2.4.gtf"

library(Rsubread)
fc <- featureCounts(bam, 
                    annot.ext = gtf, 
                    isGTFAnnotationFile = TRUE, 
                    isPairedEnd = TRUE, 
                    nthreads=8)

library(edgeR)
library(limma)
x <- DGEList(counts=fc$counts, genes=fc$annotation, group=factor(c(rep("H1975", 2), rep("HCC", 2))))
dim(x)
colnames(x) <- libs
x$samples
filter <- filterByExpr(x)
table(filter)
x <- x[filter,]
x <- calcNormFactors(x)
x$samples
cpm <- cpm(x, log=TRUE)
plotMDS(cpm, col=c("red", "blue")[as.numeric(x$samples$group)], pch = 1)
legend("topright", c("H1975", "HCC827"), text.col=c("red", "blue"), pch = 1, col = c("red", "blue"))
design <- model.matrix((~x$samples$group))
v <- voom(x, design=design, plot=TRUE, save.plot=TRUE)
fit <- lmFit(v, design = design)
efit <- eBayes(fit)
dt <- decideTests(efit)
summary(dt)

tt <- topTable(efit, coef=2, number=Inf)
anno <- read.table("../annotation/sequins/rnasequin_genes_2.4.tsv", header = TRUE, stringsAsFactors = FALSE)
anno$logFC <- log(anno$MIX_B / anno$MIX_A)
tt$logFC_expected <- anno$logFC[match(rownames(tt), anno$NAME)]
tt$LENGTH <- anno$LENGTH[match(rownames(tt), anno$NAME)]
cor(tt$logFC, tt$logFC_expected)
lfc.lm <-lm(tt$logFC ~ tt$logFC_expected)
summary(lfc.lm)

table(prediction=(tt$adj.P.Val<0.05), trueLabel=tt$logFC_expected!=0)
# explore false positive
tt[tt$adj.P.Val<0.05 & tt$logFC_expected==0,]
# explore false negative
tt[tt$adj.P.Val>0.05 & tt$logFC_expected!=0,]

write.table(tt, file = "topTable.txt", sep = "\t")