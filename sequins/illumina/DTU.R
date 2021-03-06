datadir="/stornext/General/data/user_managed/grpu_mritchie_1/long_RNA_benchmark/Illumina_April19/results/salmon"
samples <- c("1-H1975-1_S1","2-H1975-2_S2", "6-HCC-1_S6", "7-HCC-2_S7")
quant <- file.path(datadir, samples, "quant.sf")
library(tximport)
txi <- tximport(quant, type="salmon", txOut=TRUE, countsFromAbundance = "scaledTPM")
library(GenomicFeatures)
gtf <- "../flames/results/isoform_annotated.filtered.gff3"
txdb <- makeTxDbFromGFF(gtf)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID") 
tab <- table(txdf$GENEID) 
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

library(DRIMSeq)
samples <- data.frame(
  sample_id = make.names(samples),
  group = c("A", "A", "B", "B"),
  stringsAsFactors = FALSE
)
counts <- as.data.frame(txi$counts, stringAsFactors = FALSE)
counts <- counts[match(txdf$TXNAME, rownames(counts)),]
counts <- cbind(txdf$GENEID, txdf$TXNAME, counts)
colnames(counts) <- c("gene_id", "feature_id", samples$sample_id)
d <- dmDSdata(counts = counts, samples = samples)
d

anno <- read.delim("../annotations/rnasequin_isoforms_2.4.tsv", sep = "\t", stringsAsFactors = FALSE)
m <- match(anno$NAME, as.character(counts$feature_id))
annoData <- cbind(anno, counts[m, 3:6])
annoData <- data.table::melt(annoData, id.vars = 1:4)
annoData$exp <-annoData$MIX_A
annoData$exp[annoData$variable %in% c("6-HCC-1_S6", "7-HCC-2_S7")] <- 
  annoData$MIX_B[annoData$variable %in% c("6-HCC-1_S6", "7-HCC-2_S7")]
colnames(annoData)[5] <- "sample"
cor(annoData$exp, annoData$value, use = "complete.obs")

library(ggplot2)
library(scales)
# Supp Fig S7
pdf("plots/countsVsAnno.pdf", height = 5, width = 8)
ggplot(annoData, aes(x=exp, y=value, colour=sample)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(trans = "log10", labels = comma) +
  scale_y_continuous(trans = "log10", labels=comma) +
  labs(x = "Expected abundance", y = "TPM") +
  scale_colour_manual(values = c("#869700", "#D12300", "#071535", "#54B5E1"))+
  theme(text = element_text(size = 20)) +
  annotate(geom="text", x=100, y=3000000, label=paste("Pearson's r=", round(cor(annoData$exp, annoData$value, use = "complete.obs"), 2)), size=8) 
dev.off()

d <- dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_expr = 2, min_gene_expr = 10, min_feature_expr = 10, run_gene_twice = T)
plotData(d)
design <- model.matrix(~group, data=DRIMSeq::samples(d))
colnames(design) <- sub("group", "", colnames(design))
set.seed(1904)
d <- suppressWarnings(dmPrecision(d, design, BPPARAM = BiocParallel::MulticoreParam(8)))
plotPrecision(d)
d <- suppressWarnings(dmFit(d, design = design, verbose=1, BPPARAM = BiocParallel::MulticoreParam(16)))
head(proportions(d))
d <- dmTest(d, coef = "B")
head(results(d))
head(results(d, level = "feature"))

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

library(edgeR)
library(limma)
y <- DGEList(counts = as.matrix(DRIMSeq::counts(d)[,c(-1, -2)]), samples = DRIMSeq::samples(d))
rownames(y) <- counts(d)$feature_id
y <- calcNormFactors(y)
y$samples
cpm <- cpm(y, log=T)
plotMDS(cpm, col=as.numeric(as.factor(y$samples$group)))
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

pScreen <- ts$P.Value
pScreen2 <- ts.f$P.Value
names(pScreen) <- ts$GeneID
names(pScreen2) <- ts.f$GeneID
pConfirmation <- matrix(ts.tx$P.Value, ncol = 1)
rownames(pConfirmation) <- ts.tx$ExonID
tx2gene <- ts.tx[,1:2]
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

proportion <- readRDS("../proportion.RDS")
a <- match(unique(drim.padj$geneID), drim.padj$geneID)
b <- match(unique(ds.padj$geneID), ds.padj$geneID)
c <- match(unique(ds.padj.f$geneID), ds.padj.f$geneID)
results.all <- data.frame(
  id = c(
    res$gene_id,
    res.txp$feature_id,
    unique(as.character(drim.padj$geneID)),
    as.character(drim.padj$txID),
    ts$GeneID,
    ts.f$GeneID,
    ts.tx$ExonID,
    unique(as.character(ds.padj$geneID)),
    as.character(ds.padj$txID)
  ),
  FDR = c(
    res$adj_pvalue,
    res.txp$adj_pvalue,
    drim.padj$gene[a],
    drim.padj$transcript,
    ts$FDR,
    ts.f$FDR,
    ts.tx$FDR,
    ds.padj$gene[b],
    ds.padj$transcript
  ),
  DTU = c(
    proportion$isChanges[match(res$gene_id, proportion$geneID)],
    proportion$isChanges[match(res.txp$feature_id, proportion$txID)],
    proportion$isChanges[match(drim.padj$txID, proportion$txID)][a],
    proportion$isChanges[match(drim.padj$txID, proportion$txID)],
    proportion$isChanges[match(ts$GeneID, proportion$geneID)], 
    proportion$isChanges[match(ts.f$GeneID, proportion$geneID)],
    proportion$isChanges[match(ts.tx$ExonID, proportion$txID)], 
    proportion$isChanges[match(ds.padj$geneID, proportion$geneID)][b], 
    proportion$isChanges[match(ds.padj$txID, proportion$txID)]
  ),
  test = c(
    rep("DRIMSeq_gene", nrow(res)),
    rep("DRIMSeq_feature", nrow(res.txp)),
    rep("DRIMSeq_stageR_gene", length(unique(drim.padj$geneID))),
    rep("DRIMSeq_stageR_transcript", nrow(drim.padj)),
    rep("diffSplice_simes", nrow(ts)),
    rep("diffSplice_F", nrow(ts.f)),
    rep("diffSplice_t", nrow(ts.tx)),
    rep("diffSplice_stageR_gene", length(unique(ds.padj$geneID))),
    rep("diffSplice_stageR_transcript", nrow(ds.padj))
  ),
  level = c(
    rep("Gene level", nrow(res)),
    rep("Transcript level", nrow(res.txp)),
    rep("Gene level", length(unique(drim.padj$geneID))),
    rep("Transcript level", nrow(drim.padj)),
    rep("Gene level", nrow(ts)),
    rep("Gene level", nrow(ts.f)),
    rep("Transcript level", nrow(ts.tx)),
    rep("Gene level", length(unique(ds.padj$geneID))),
    rep("Transcript level", nrow(ds.padj))
  )
)

results.all$test <-
  factor(
    results.all$test,
    levels = c(
      "DRIMSeq_gene",
      "DRIMSeq_stageR_gene",
      "diffSplice_simes",
      "diffSplice_F",
      "diffSplice_stageR_gene",
      "DRIMSeq_feature",
      "DRIMSeq_stageR_transcript",
      "diffSplice_t",
      "diffSplice_stageR_transcript"
    )
  )

# manual colour scheme
col_impression <- c("#D96A70", "#C68C81", "#303C24", "#476937", "#9FC675", "#8E3055", "#D5A2CB", "#3C5880", "#708FA6")

# Supp Fig S8A
p1 <- ggplot(na.omit(results.all), aes(x = DTU, y = FDR, fill=test, colour=test)) +
  geom_boxplot(alpha=.6, outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.height = .02), cex=.7, alpha = .7)+
  facet_grid(.~level) +
  theme_bw() +
  theme(text = element_text(size=20), legend.text = element_text(size=14), legend.position = "none")+
  geom_hline(yintercept = 0.05, linetype="dashed") +
  labs(y = "Adjusted p-value", x="is DTU") +
  scale_fill_manual(values = col_impression) +
  scale_colour_manual(values = col_impression)

results.all$sig <- results.all$FDR < 0.05
FDR <- sapply(unique(results.all$test), function(x){
  results = na.omit(results.all[results.all$test == x, ])
  fd = sum(results$sig & !(results$DTU))
  dt = sum(results$sig)
  return(c(as.character(x), fd/dt))
}, simplify = TRUE)
FDR <- as.data.frame(t(FDR))
FDR$V2 <- as.numeric(as.character(FDR$V2))
colnames(FDR) <- c("test", "FDR")
FDR$test <- factor(FDR$test, levels = c("DRIMSeq_gene", "DRIMSeq_stageR_gene","diffSplice_simes", "diffSplice_F", "diffSplice_stageR_gene", "DRIMSeq_feature", "DRIMSeq_stageR_transcript",
                               "diffSplice_t", "diffSplice_stageR_transcript"))
FDR

TPR <- sapply(unique(results.all$test), function(x){
  results = na.omit(results.all[results.all$test == x, ])
  td = sum(results$sig & results$DTU)
  pc = sum(results$DTU)
  return(c(as.character(x), td/pc))
}, simplify = TRUE)
TPR <- as.data.frame(t(TPR))
TPR$V2 <- as.numeric(as.character(TPR$V2))
colnames(TPR) <- c("test", "TPR")
TPR$test <- factor(TPR$test, levels = c("DRIMSeq_gene", "DRIMSeq_stageR_gene","diffSplice_simes", "diffSplice_F", "diffSplice_stageR_gene", "DRIMSeq_feature", "DRIMSeq_stageR_transcript",
                               "diffSplice_t", "diffSplice_stageR_transcript" ))
TPR


# Supp Fig S8C
p3 <- ggplot(TPR, aes(x=test, y=TPR, fill=test)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs (x = "Test", y = "True positive rate") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), text = element_text(size=16),
        legend.position = "none") +
  scale_fill_manual(values = col_impression)
# Supp fig S8B
p2 <- ggplot(FDR, aes(x=test, y=FDR, fill=test)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs (x = "Test", y = "False discovery rate") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), text = element_text(size=16),
        legend.position = "none") +
  scale_fill_manual(values = col_impression)

library(cowplot)
pdf("plots/DTU_short.pdf", height = 10, width = 8)
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 24, nrow=3)
dev.off()

FP <- na.omit(results.all[results.all$DTU==FALSE & results.all$FDR<0.05,])
FN <- na.omit(results.all[results.all$DTU==TRUE & results.all$FDR>0.05,])
FP
FN

# gene level diffSplice result (simes)
write.table(ts, "topSplice_simes.tsv", sep = "\t")
# drimseq+stager results
write.table(drim.padj, "drimseq_stager.tsv", sep = "\t")
# drimseq gene
write.table(res, "drimseq_gene.tsv", sep = "\t")
# drimseq proportion
write.table(proportions(d), "drimseq_proportion.tsv", sep="\t")