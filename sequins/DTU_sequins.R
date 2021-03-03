library(limma)
library(edgeR)
library(data.table)
library(ggplot2)
library(scales)
library(DRIMSeq)
library(stageR)

data <- read.csv("./flames/results/transcript_count.csv", stringsAsFactors = FALSE)
anno <- read.delim("./annotations/rnasequin_isoforms_2.4.tsv", sep = "\t", stringsAsFactors = FALSE)
geneanno <- read.delim("./annotations/rnasequin_genes_2.4.tsv", sep = "\t", stringsAsFactors = FALSE)
# calculate annotate proportion
proportion <- data.frame(txID = anno$NAME)
tx2gene <- limma::strsplit2(proportion$txID, "_")
tx2gene <-  paste(tx2gene[,1], tx2gene[,2], sep="_")
proportion$geneID <- tx2gene
m <- match(proportion$geneID, geneanno$NAME)
proportion$MIX_A <- anno$MIX_A / geneanno$MIX_A[m]
proportion$MIX_B <- anno$MIX_B / geneanno$MIX_B[m]
proportion$FC <- proportion$MIX_A / proportion$MIX_B
proportion$isChanges <- FALSE
proportion$isChanges[abs(proportion$FC-1) > 0.001] <- TRUE

cpm <- cpm(data[3:6])
m <- match(anno$NAME, data$transcript_id)
annoData <- cbind(anno, cpm[m,])
colnames(annoData)[5:8] <- paste0("mix", c("B1", "A2", "B2", "A1"))
annoData <- melt(annoData, id.vars = 1:4)
annoData$exp <-annoData$MIX_A
annoData$exp[annoData$variable %in% c("mixB1", "mixB2")] <- annoData$MIX_B[annoData$variable %in% c("mixB1", "mixB2")]
colnames(annoData)[5] <- "sample"
annoData$sample <- factor(annoData$sample, levels = c("mixA1", "mixA2", "mixB1", "mixB2"))
cor(annoData$exp, annoData$value, use = "complete.obs")
# Fig 3D
pdf("plots/countsVsAnno.pdf", height = 5, width = 8)
ggplot(annoData, aes(x=exp, y=value, colour=sample)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(trans = "log10", labels = comma) +
  scale_y_continuous(trans = "log10", labels=comma, limits=c(0.5, 180000)) +
  labs(x = "Expected abundance", y = "CPM") +
  scale_colour_manual(values = c("#869700", "#D12300", "#071535", "#54B5E1"))+
  theme(text = element_text(size = 20)) +
  annotate(geom="text", x=300, y=160000, label=paste0("Pearson's r=", round(cor(annoData$exp, annoData$value, use = "complete.obs"), 2)), size=8) 
dev.off()

#----
# DTU using DRIMSeq
counts <- data[,c(2, 1, 3, 4, 5, 6)]
colnames(counts)[2] <- "feature_id"
samples <- data.frame(
  sample_id = colnames(counts)[3:6],
  group = c("B", "A", "B", "A")
)
d <- dmDSdata(counts=counts, samples = samples)
d <- dmFilter(d, min_samps_gene_expr = 4, min_samps_feature_expr = 2, min_gene_expr = 10, min_feature_expr = 10, run_gene_twice = T)

design <- model.matrix(~group, data=DRIMSeq::samples(d))
colnames(design) <- sub("group", "", colnames(design))
set.seed(1904)
d <- suppressWarnings(dmPrecision(d, design, BPPARAM = BiocParallel::MulticoreParam(16)))
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

# stageR
## Assign gene-level pvalues to the screening stage
pScreen <- res$pvalue
# strp <- function(x) substr(x,1,18)
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
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=FALSE)
})
head(drim.padj)

#----
# DTU using diffSplice
y <- DGEList(counts = as.matrix(DRIMSeq::counts(d)[,c(-1, -2)]), samples = DRIMSeq::samples(d))
rownames(y) <- counts(d)$feature_id
y <- calcNormFactors(y)
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

# stageR
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
  ds.padj <- getAdjustedPValues(stageRObj2, order=FALSE,
                                  onlySignificantGenes=FALSE)
  ds.padj.f <- getAdjustedPValues(stageRObj3, order=FALSE,
                                  onlySignificantGenes=FALSE)
})
head(ds.padj)
head(ds.padj.f)

#----
# Explore DTU results
# drimseq feature level
res.tx <- results(d, level = "feature")
res.tx$isChanges <- proportion$isChanges[match(res.tx$feature_id, proportion$txID)]
# drimseq gene level
res.gene <- results(d, level = "gene")
res.gene$isChanges <- proportion$isChanges[match(res.gene$gene_id, proportion$geneID)]
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=TRUE,
                                  onlySignificantGenes=FALSE)
})

drim.padj$isChanges <- proportion$isChanges[match(drim.padj$txID, proportion$txID)]
a <- match(unique(drim.padj$geneID), drim.padj$geneID)
b <- match(unique(ds.padj$geneID), ds.padj$geneID)
c <- match(unique(ds.padj.f$geneID), ds.padj.f$geneID)
results.all <- data.frame(
  id = c(
    res.gene$gene_id,
    res.tx$feature_id,
    unique(as.character(drim.padj$geneID)),
    as.character(drim.padj$txID),
    ts$GeneID,
    ts.f$GeneID,
    ts.tx$ExonID,
    unique(as.character(ds.padj$geneID)),
    as.character(ds.padj$txID)
  ),
  FDR = c(
    res.gene$adj_pvalue,
    res.tx$adj_pvalue,
    drim.padj$gene[a],
    drim.padj$transcript,
    ts$FDR,
    ts.f$FDR,
    ts.tx$FDR,
    ds.padj$gene[b],
    ds.padj$transcript
  ),
  DTU = c(
    res.gene$isChanges,
    res.tx$isChanges,
    drim.padj$isChanges[a],
    drim.padj$isChanges,
    proportion$isChanges[match(ts$GeneID, proportion$geneID)],
    proportion$isChanges[match(ts.f$GeneID, proportion$geneID)],
    proportion$isChanges[match(ts.tx$ExonID, proportion$txID)],
    proportion$isChanges[match(ds.padj$geneID, proportion$geneID)][b],
    proportion$isChanges[match(ds.padj$txID, proportion$txID)]
  ),
  test = c(
    rep("DRIMSeq_gene", nrow(res.gene)),
    rep("DRIMSeq_feature", nrow(res.tx)),
    rep("DRIMSeq_stageR_gene", length(unique(drim.padj$geneID))),
    rep("DRIMSeq_stageR_transcript", nrow(drim.padj)),
    rep("diffSplice_simes", nrow(ts)),
    rep("diffSplice_F", nrow(ts.f)),
    rep("diffSplice_t", nrow(ts.tx)),
    rep("diffSplice_stageR_gene", length(unique(ds.padj$geneID))),
    rep("diffSplice_stageR_transcript", nrow(ds.padj))
  ),
  level = c(
    rep("Gene level", nrow(res.gene)),
    rep("Transcript level", nrow(res.tx)),
    rep("Gene level", length(unique(drim.padj$geneID))),
    rep("Transcript level", nrow(drim.padj)),
    rep("Gene level", nrow(ts)),
    rep("Gene level", nrow(ts.f)),
    rep("Transcript level", nrow(ts.tx)),
    rep("Gene level", length(unique(ds.padj$geneID))),
    rep("Transcript level", nrow(ds.padj))
  )
)
results.all$test <- factor(results.all$test, levels = c("DRIMSeq_gene", "DRIMSeq_stageR_gene","diffSplice_simes", "diffSplice_F", "diffSplice_stageR_gene", "DRIMSeq_feature", "DRIMSeq_stageR_transcript",
                               "diffSplice_t", "diffSplice_stageR_transcript"))
col_impression <- c("#D96A70", "#C68C81", "#303C24", "#476937", "#9FC675", "#8E3055", "#D5A2CB", "#3C5880", "#708FA6")
# Supp fig S6
pdf("plots/DTU_FDR.pdf", width = 17, height = 5)
ggplot(na.omit(results.all), aes(x = DTU, y = FDR, fill=test, colour=test)) +
  geom_boxplot(alpha=.6, outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.height = .02), cex=.7, alpha = .7)+
  facet_grid(.~level) +
  theme_bw() +
  theme(text = element_text(size=20), legend.text = element_text(size=14))+
  geom_hline(yintercept = 0.05, linetype="dashed") +
  labs(y = "Adjusted p-value", x="is DTU") +
  scale_fill_manual(values = col_impression) +
  scale_colour_manual(values = col_impression)
dev.off()

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

# Fig 3E and 3F
pdf("plots/DTU_TPR.pdf", height=5, width=8)
ggplot(TPR, aes(x=test, y=TPR, fill=test)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs (x = "Test", y = "True positive rate") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), text = element_text(size=16),
        legend.position = "none") +
  scale_fill_manual(values = col_impression)
dev.off()

pdf("plots/DTU_FDR0.pdf", height=5, width=8)
ggplot(FDR, aes(x=test, y=FDR, fill=test)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs (x = "Test", y = "False discovery rate") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), text = element_text(size=16),
        legend.position = "none") +
  scale_fill_manual(values = col_impression)
dev.off()

FP <- na.omit(results.all[results.all$DTU==FALSE & results.all$FDR<0.05,])
FN <- na.omit(results.all[results.all$DTU==TRUE & results.all$FDR>0.05,])

#----
# Compare to short-read
ts.short <- read.table("./illumina/topSplice_simes.tsv", stringsAsFactors = FALSE, sep = "\t", header=TRUE)
calcTopN <- function(tl, ts, n){
  num <- sapply(1:n, function(x){
    length(intersect(tl[1:x, "GeneID"], ts[1:x, "GeneID"]))
  }, simplify = TRUE)
  intTopN <- data.frame(
    number = 1:n,
    intersect = num
  )
  return(intTopN)
}
top47 <- calcTopN(ts, ts.short, 47)

res.short <- read.table("./illumina/drimseq_gene.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
calcTopN_drim <- function(res.long, res.short, n){
  # sort results
  res.long <- res.long[order(res.long$adj_pvalue), ]
  res.short <- res.short[order(res.short$adj_pvalue), ]
  num <- sapply(1:n, function(x){
    length(intersect(res.long[1:x, "gene_id"], res.short[1:x, "gene_id"]))
  }, simplify = TRUE)
  intTopN <- data.frame(
    number = 1:n,
    intersect = num
  )
  return(intTopN)
}
top47.drim <- calcTopN_drim(res, res.short, 47)

top47.all <- data.frame(
  number = c(top47$number, top47.drim$number),
  intersect = c(top47$intersect, top47.drim$intersect),
  method = rep(c("diffSplice", "DRIMSeq"), c(47, 47))
)
# Supp fig S9
pdf("plots/topN.pdf", height = 5, width = 8)
ggplot(top47.all, aes(x=number, y=intersect, colour=method)) + geom_line() +
  geom_abline(colour="grey") +
  theme_bw() +
  theme(text = element_text(size=20), legend.text = element_text(size=14)) +
  labs(x = "Top n genes", y = "Number of intersect") +
  scale_color_brewer(palette="Set1")
dev.off()



proportions.short <- read.table("illumina/drimseq_proportion.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
# only compare wt proportion
m <- match(proportions.short$feature_id, proportions(d)$feature_id)
cor(proportions(d)[m, 3], proportions.short$X6.HCC.1_S6, use = "complete.obs")
plot(proportions(d)[m, 3], proportions.short$X6.HCC.1_S6)