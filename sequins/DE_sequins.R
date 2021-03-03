library(Rsubread)
library(limma)
library(edgeR)
library(ggplot2)
library(scales)

bams <- paste0("./talon/labeled/barcode0", 1:4, "_clean_labeled.sam")
gff <- "/wehisan/home/allstaff/d/dong.x/annotation/sequins/rnasequin_annotation_2.4.gtf"
fc <- featureCounts(bams, annot.ext=gff, isGTFAnnotationFile = TRUE, isLongRead = TRUE, primaryOnly=TRUE)
fc$stat

x <- DGEList(counts = fc$counts)
x$samples$group <- c("A", "A", "B", "B")
filter <- filterByExpr(x)
summary(filter)
x <- x[filter, keep.lib.sizes=FALSE]
x <- calcNormFactors(x)
lcpm <- cpm(x, log=TRUE)
plotMDS(lcpm, labels = c("mixA1", "mixA2", "mixB1", "mixB2"), col = c("#071535", "#071535", "#D12300", "#D12300"))

design <- model.matrix(~x$samples$group)
colnames(design) <- sub("x\\$samples\\$group", "", colnames(design))
colnames(design) <- sub("x\\$samples\\$", "", colnames(design))
# Supp fig S3
pdf("plots/voom_sequins.pdf", height = 4)
v <- voomWithQualityWeights(x, design=design, plot=TRUE)
dev.off()

fit <- lmFit(v, design = design)
efit <- eBayes(fit)
dt <- decideTests(efit)
summary(dt)

plotMD(efit, status=dt[,2], column = 2, values=c(1,-1), hl.col=c("red", "blue"),
       main = "mix B vs A")

tt <- topTable(efit, coef=2, number=Inf)
anno <- read.table("./annotations/rnasequin_genes_2.4.tsv", header = TRUE, stringsAsFactors = FALSE)
anno$logFC <- log(anno$MIX_B / anno$MIX_A)
tt$logFC_expected <- anno$logFC[match(rownames(tt), anno$NAME)]
tt$LENGTH <- anno$LENGTH[match(rownames(tt), anno$NAME)]
cor(tt$logFC, tt$logFC_expected)
lfc.lm <-lm(tt$logFC ~ tt$logFC_expected)
summary(lfc.lm)
# Fig 2A
pdf("plots/logfc.pdf", height = 5, width = 8)
ggplot(tt, aes(x = logFC_expected, y = logFC)) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 20), 
        legend.text = element_text(size=16)) +
  labs(x = "Expected logFC", y="Estimated logFC") +
  scale_x_continuous(breaks = seq(-4, 4, 2), limits = c(-4.1, 4.1)) +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4.4, 5.6)) +
  geom_smooth(method="lm") +
  annotate(geom = "text", x = 3, y = 5.3, label = paste("Adj R2 = ",round(summary(lfc.lm)$adj.r.squared, 2)), size = 7)
dev.off()
# Supp Fig S2
cor(tt$AveExpr, tt$LENGTH)
pdf("plots/lengthAveExpSequins.pdf")
ggplot(tt, aes(x = LENGTH, y = AveExpr)) +
  geom_point() +
  theme_bw() +
  labs(x = "Gene length", y = "Average expression (log-CPM)") +
  scale_x_continuous(trans = "log10") +
  theme(text = element_text(size = 16))
dev.off()

table(prediction=(tt$adj.P.Val<0.05), trueLabel=tt$logFC_expected!=0)
tt$sig <- tt$adj.P.Val<0.05
tt$DE <- abs(tt$logFC_expected) > 0.0001
FDR = sum(tt$sig &!(tt$DE)) / sum(tt$sig)
TPR = sum(tt$DE & tt$sig) / sum(tt$DE)

# Compare against short-read results
tt.short <-read.delim("./illumina/topTable.txt", sep = "\t", stringsAsFactors = FALSE)
tt$t.short <- tt.short$t[match(rownames(tt), rownames(tt.short))]
                         
t.lm <- lm(tt$t.short ~ tt$t)
summary(t.lm)
# Fig 2B
pdf("plots/t.pdf", height = 5, width = 8)
ggplot(tt, aes(x = t, y = t.short))+
  geom_point() + 
  labs(x = "Long read t-statistic", y = "Short read t-statistic") +
  geom_smooth(method='lm', formula= y~x)+
  theme_bw() +
  theme(text = element_text(size = 20), 
        legend.text = element_text(size=16)) +
  annotate(geom = "text", x = 13, y = 90, label = paste("Adj R2 = ",round(summary(t.lm)$adj.r.squared, 2)), size = 7)
dev.off()

# Supp fig S1, read number come from sequencing summary file.
read.num <- data.frame(
  category = factor(c("raw reads", "demultiplexed reads", "gene counts"), levels = c("raw reads", "demultiplexed reads", "gene counts")),
  Number = c(7349818, 7053814, 5684786)
)
library(scales)
pdf("readnumsequins.pdf")
ggplot(read.num, aes(x=category, y=Number)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  labs(x="")+
  geom_text(aes(label=Number), position = position_stack(vjust = 0.5), colour="white")
dev.off()

x$samples$sample <- c("mixA1", "mixA2", "mixB1", "mixB2")
ggplot(x$samples, aes(x=sample, y=lib.size, fill=group)) +
  geom_bar(stat="identity")+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(y = "Library size", x = "Sample")