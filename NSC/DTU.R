library(DRIMSeq)
library(limma)
library(edgeR)

scBarcode <- read.delim("./flames/results/pseudo_barcode_annotation.csv", sep = ",", stringsAsFactors = FALSE)
scBarcode$barcode <- sub(".fq.gz", "", scBarcode$file_name)
scBarcode <- scBarcode[order(scBarcode$barcode), ]
scBarcode$group <- c("WT", "Smchd1-null", "Smchd1-null", "WT", "WT",  "WT", "WT", "Smchd1-null")
scBarcode

counts <- read.csv("./flames/results/transcript_count.csv")
m <- match(colnames(counts)[3:ncol(counts)], scBarcode$pseudo_barcode)
colnames(counts)[3:ncol(counts)] <- scBarcode$barcode[m]
colnames(counts)[1] <- "feature_id"
# combine the wt sample that was sequenced twice
counts$barcode07 <- counts$barcode07 + counts$barcode13
BC13 <- which(colnames(counts)=="barcode13")
counts <- counts[,-BC13]
colnames(counts)[3:ncol(counts)] <- make.names(colnames(counts)[3:ncol(counts)])
samples <- data.frame(sample_id = colnames(counts)[3:ncol(counts)], group=scBarcode$group[match(colnames(counts)[3:ncol(counts)], scBarcode$barcode)])

geneCount <- aggregate(counts[,c(-1,-2)], by=list(counts$gene_id), FUN=sum)
rownames(geneCount) <- geneCount[,1]
geneCount <- geneCount[,-1]
fc <- readRDS("counts.RDS")
fcCount <- fc$counts
colnames(fcCount) <- sub(".sorted.bam", "", colnames(fcCount))
fcCount[,"barcode07"] <- fcCount[,"barcode07"] + fcCount[,"barcode13"]
fcCount <- fcCount[, -5]
samplenames = paste("sample", c(5, 7, 3, 6, 2, 4, 1))

## ----DRIMSeq--------------------------------------------------------
d <- dmDSdata(counts = counts, samples = samples)
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, min_gene_expr = 10, min_feature_expr = 10, run_gene_twice = T)
plotData(d)

design <- model.matrix(~group, data=DRIMSeq::samples(d))
colnames(design) <- sub("group", "", colnames(design))
set.seed(1904)
d <- suppressWarnings(dmPrecision(d, design, BPPARAM = BiocParallel::MulticoreParam(16)))

head(mean_expression(d))
common_precision(d)
head(genewise_precision(d))

d <- suppressWarnings(dmFit(d, design = design, verbose=1, BPPARAM = BiocParallel::MulticoreParam(16)))
head(proportions(d))

d <- dmTest(d, coef = "WT")
head(results(d))
head(results(d, level = "feature"))

res <- results(d)
res.txp <- results(d, level = "feature")
## turn NA p value into 1
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)
write.table(res, "DTUres.tsv", sep = "\t", row.names = FALSE)

library(RColorBrewer)
library(ggplot2)
col <- brewer.pal(3, "Set2")

library(stageR)
## Assign gene-level pvalues to the screening stage
pScreen <- res$pvalue
names(pScreen) <- results(d)$gene_id
## Assign transcript-level pvalues to the confirmation stage
pConfirmation <- matrix(res.txp$pvalue, ncol = 1)
rownames(pConfirmation) <- res.txp$feature_id
## Create the gene-transcript mapping
tx2gene <- res.txp[,c("feature_id", "gene_id")]
## Create the stageRTx object and perform the stage-wise analysis
stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
alpha = 0.25)

suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=TRUE,
                                  onlySignificantGenes=FALSE)
})
head(drim.padj)
write.table(drim.padj, file = "DTUpadj.txt", sep = "\t", row.names = FALSE)

pdf("plots/DTUPropPisd.pdf", height = 6, width =9)
plotProportions(d, gene_id = "ENSMUSG00000023452.19", group_variable = "group", 
               group_colors = col[1:2], plot_type = "boxplot1") + 
  ggtitle("Pisd") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5), 
        text = element_text(size=10), 
        legend.text=element_text(size=10)) 
dev.off()

## ----diffSplice----------------------------------------------------------
library(edgeR)
library(limma)
y <- DGEList(counts = as.matrix(DRIMSeq::counts(d)[,c(-1, -2)]), samples = DRIMSeq::samples(d))
rownames(y) <- counts(d)$feature_id
y <- calcNormFactors(y)
y$samples
cpm <- cpm(y, log=T)
library(Mus.musculus)
geneid <- counts(d)$gene_id
geneid <- substr(geneid, 1, 18)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM","ENTREZID"), 
                keytype="ENSEMBL")
m <- match(geneid, genes$ENSEMBL)
genes <- genes[m,]
head(genes)
y$genes <- genes

v <- voom(y, design, plot=T)
fit <- lmFit(v,design)
efit <- eBayes(fit)
geneid <- counts(d)$gene_id
featureid <- counts(d)$feature_id

ex <- diffSplice(efit, geneid = geneid, exonid=featureid)
ts <- topSplice(ex,  number=Inf)
ts.tx <- topSplice(ex, test = "t", number = Inf)
# plotSplice(ex, geneid = "Pisd", genecolname="SYMBOL", FDR = 0.25)

library(stageR)

pScreen <- ts$P.Value
names(pScreen) <- ts$GeneID
pConfirmation <- matrix(ts.tx$P.Value, ncol = 1)
rownames(pConfirmation) <- ts.tx$ExonID

tx2gene <- ts.tx[, c(6, 5)]

stageRObj2 <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj2 <- stageWiseAdjustment(object = stageRObj2, method = "dtu",
alpha = 0.25)
suppressWarnings({
  ds.padj <- getAdjustedPValues(stageRObj2, order=TRUE,
                                  onlySignificantGenes=FALSE)
})
head(ds.padj)

# plot Pisd gene model
library(Gviz)
library(GenomicFeatures)

gff <- "./flames/results/isoform_annotated.filtered.gff3"
anno <- makeTxDbFromGFF(gff, organism = "Mus musculus")

tr_Pisd <- GeneRegionTrack(anno, chromosome = "chr5", 
                           start = 32736314, end = 32785626)
gtrack <- GenomeAxisTrack()
pdf("plots/genePisd.pdf", width = 10, height = 5)
plotTracks(list(gtrack, tr_Pisd), lwd = 1, lwd.grid = 1)
dev.off()

save.image("DTU.RData")