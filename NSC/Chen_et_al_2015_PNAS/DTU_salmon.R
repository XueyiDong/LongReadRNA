# Load Salmon quant TPM data
datadir="/wehisan/general/academic/seq_data/kelan/AGRF_CAGRF6583/salmon"
samples <- c("113_C26VPACXX_ATCACG", "114_C26VPACXX_CGATGT", "11_C26VPACXX_GATCAG", "12_C26VPACXX_TAGCTT", "13_C26VPACXX_GGCTAC", "14_C26VPACXX_CTTGTA")
quant <- file.path(datadir, samples, "quant.sf")
library(tximport)
txi <- tximport(quant, type="salmon", txOut=TRUE, countsFromAbundance = "scaledTPM")
library(GenomicFeatures)
gtf <- "../flames/results/isoform_annotated.filtered.gff3"
txdb <- makeTxDbFromGFF(gtf)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID") 
tab <- table(txdf$GENEID) 
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

## -------------------------------------------------------
library(DRIMSeq)
samples <- data.frame(
  sample_id = samples,
  group = c("WT", "WT", "Smchd1-null", "Smchd1-null", "Smchd1-null", "WT"),
  stringsAsFactors = FALSE
)
samples$group <- as.factor(samples$group)
counts <- as.data.frame(txi$counts, stringAsFactors = FALSE)
counts <- counts[match(txdf$TXNAME, rownames(counts)),]
counts <- cbind(txdf$GENEID, txdf$TXNAME, counts)
colnames(counts)[3:ncol(counts)] <- make.names(colnames(counts)[3:ncol(counts)])
colnames(counts) <- c("gene_id", "feature_id", samples$sample_id)
d <- dmDSdata(counts = counts, samples = samples)
d <- dmFilter(d, min_samps_gene_expr = 6, min_samps_feature_expr = 3, min_gene_expr = 10, min_feature_expr = 10, run_gene_twice = T)

design <- model.matrix(~group, data=DRIMSeq::samples(d))
colnames(design) <- sub("group", "", colnames(design))
set.seed(1904)
d <- suppressWarnings(dmPrecision(d, design, BPPARAM = BiocParallel::MulticoreParam(8)))
d <- suppressWarnings(dmFit(d, design = design, verbose=1, BPPARAM = BiocParallel::MulticoreParam(16)))
head(proportions(d))

d <- dmTest(d, coef = "WT")
head(results(d))
head(results(d, level = "feature"))

library(RColorBrewer)
library(ggplot2)
col <- brewer.pal(3, "Set2")
pdf("shortDTUPropPisd.pdf", height = 6, width =9)
plotProportions(d, gene_id = "ENSMUSG00000023452.19", group_variable = "group", 
               group_colors = col[1:2], plot_type = "boxplot1") + 
  ggtitle("Pisd (short read)") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5), 
        text = element_text(size=10), 
        legend.text=element_text(size=10)) 
dev.off()

res <- results(d)
res.txp <- results(d, level = "feature")
no.na <- function(x) ifelse(is.na(x), 1, x) 
res$pvalue <- no.na(res$pvalue) 
res.txp$pvalue <- no.na(res.txp$pvalue)
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
alpha = 0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=TRUE,
                                  onlySignificantGenes=FALSE)
})
head(drim.padj)

## -------------------------------------------------------
library(edgeR)
library(limma)
y <- DGEList(counts = as.matrix(DRIMSeq::counts(d)[,c(-1, -2)]), samples = DRIMSeq::samples(d))
rownames(y) <- counts(d)$feature_id
y <- calcNormFactors(y)
y$samples
cpm <- cpm(y, log=T)
plotMDS(cpm, col=as.numeric(as.factor(y$samples$group)))

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
ts.f <- topSplice(ex, test = "F", number = Inf)
head(ts)
table(ts$FDR < 0.05)

library(stageR)
pScreen <- ts$P.Value
pScreen2 <- ts.f$P.Value
names(pScreen) <- ts$GeneID
names(pScreen2) <- ts.f$GeneID
pConfirmation <- matrix(ts.tx$P.Value, ncol = 1)
rownames(pConfirmation) <- ts.tx$ExonID
tx2gene <- ts.tx[,c("ExonID", "GeneID")]
stageRObj2 <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj2 <- stageWiseAdjustment(object = stageRObj2, method = "dtu",
alpha = 0.05)
stageRObj3 <- stageRTx(pScreen = pScreen2, pConfirmation = pConfirmation,
pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj3 <- stageWiseAdjustment(object = stageRObj3, method = "dtu",
alpha = 0.05)
suppressWarnings({
  ds.padj <- getAdjustedPValues(stageRObj2, order=TRUE,
                                  onlySignificantGenes=FALSE)
  ds.padj.f <- getAdjustedPValues(stageRObj3, order=TRUE,
                                  onlySignificantGenes=FALSE)
})
head(ds.padj)
head(ds.padj.f)
# plotSplice(ex, geneid = "Pisd", genecolname="SYMBOL")
