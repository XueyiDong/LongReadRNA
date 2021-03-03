counts <- readRDS("./counts.RDS")
library(edgeR)
library(limma)
x <- DGEList(counts = fc$counts)
colnames(x)<- c("BC07", paste0("BC", 10:13), paste0("BC", 15:17))
x$samples$group <- factor(c("WT", "Smchd1-null", "Smchd1-null", "WT", "WT",  "WT", "WT", "Smchd1-null"), levels = c("WT", "Smchd1-null"))
x$samples$run <- factor(rep(c(1, 2), c(4,4)))
x$samples

# merge from run 2 to run 1
x$counts[, 1] <- x$counts[,1] + x$counts[, 5]
x$samples$lib.size[1] <- x$samples$lib.size[1] + x$samples$lib.size[5]
x$samples$run[1] <- "NA"
# delete the sample from run 2
x <- x[,-5]
filter <- filterByExpr(x)
summary(filter)
x <- x[filter, , keep.lib.sizes = FALSE]
x <- calcNormFactors(x)
x$samples$sample <- 1:7
x$samples
lcpm <- cpm(x, log=TRUE)

library(RColorBrewer)
col <- brewer.pal(3, "Set2")
# to match color
col <- col[c(2,1,3)]
# Fig 2C
pdf("plots/mds_merged.pdf", height = 5, width = 8)
plotMDS(lcpm, col=col[x$samples$group], labels = 1:7, asp = 1,
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
legend("bottomright", c("WT", "Smchd1-null"), text.col = col[c(1,2)], cex = 1.25)
dev.off()
# Fig 1D
pdf("plots/libsize.pdf", height = 4, width = 8)
library(ggplot2)
library(scales)
ggplot(x$samples, aes(x=sample, y=lib.size, fill=group)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Set2", direction = -1) +
  scale_x_discrete(limits = x$samples$sample, labels=1:7) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  labs(y = "Library size", x = "Sample") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))
dev.off()

# Annotation
library(Mus.musculus)
geneid <- rownames(x)
geneid <- substr(geneid, 1, 18)
genes <- select(Mus.musculus, keys=geneid, 
                columns=c("SYMBOL", "TXCHROM","ENTREZID"), 
                keytype="ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),]
head(genes)
# Add gene length information, from fc
genes$length <- fc$annotation$Length[match(rownames(x), fc$annotation$GeneID)]
loglength = log(genes$length)
x$genes <- genes

design <- model.matrix(~x$samples$group)
colnames(design) <- sub("x\\$samples\\$group", "", colnames(design))
v <- voomWithQualityWeights(x, design=design, plot=TRUE, col=col[x$samples$group], save.plot=TRUE)
# Fig 2D
pdf("plots/voomTrend.pdf", height = 5, width = 8)
plot(v$voom.xy$x, v$voom.xy$y, xlab = v$voom.xy$xlab, ylab = v$voom.xy$ylab,
     pch = 16, cex = 0.25,
     cex.lab = 1.5, cex.axis = 1.5)
lines(v$voom.line, col = "red")
dev.off()
# Supp fig S4
pdf("plots/sampleWeights.pdf", height = 5, width = 8)
barplot(v$targets$sample.weights, names = 1:length(v$targets$sample.weights),
        xlab = "Sample", ylab = "Weight", col=col[x$samples$group],
        cex.lab = 1.5, cex.axis = 1.5)
abline(h = 1, col = 2, lty = 2)
legend("topright", c("WT", "Smchd1-null"), text.col = col[c(1,2)], cex = 1.25)
dev.off()

fit <- lmFit(v, design = design)
efit <- eBayes(fit)
dt <- decideTests(efit, p.value = 0.25)
summary(dt)

tt <- topTable(efit, coef=1, number=Inf)
cor(tt$length, tt$AveExpr)

library(ggplot2)
library(viridis)
# Fig 1E
pdf("plots/lengthAveExp.pdf", height = 4, width = 8)
ggplot(tt, aes(x=length, y=AveExpr)) +
  stat_binhex() +
  theme_bw() +
  labs(x="Gene length", y = "Average expression (log-CPM)", fill = "Density:\nnumber of \ngenes") +
  scale_x_continuous(trans = "log10") +
  annotate(geom="text", x=6600, y=14, label=paste0("Pearson's r=", round(cor(tt$length, tt$AveExpr), 3)), size=8) +
  scale_fill_viridis(direction = -1, option="A") +
  theme(text=element_text(size = 20)) 
dev.off()

sm <- which(efit$genes$SYMBOL=="Smchd1")
imprinted <- which(efit$genes$SYMBOL %in% c("Ndn", "Mkrn3", "Peg12"))
# Fig 2E
pdf("plots/md.pdf", height = 5, width = 8)
  plotMD(efit, status=dt[,2], column=2, values=c(1,-1), hl.col=c("red", "blue"), legend=FALSE, main=NA,
         cex.lab = 1.5, cex.axis = 1.5)
  # highlight smchd1
  points(efit$Amean[sm], efit$coefficients[sm, 2], col = "forestgreen", pch = 17, cex = 1.5)
  points(efit$Amean[imprinted], efit$coefficients[imprinted, 2], col = "goldenrod", pch = 15, cex = 1.2)
  legend("topright", legend = c("Imprinted", "Up", "NotSig", "Down", "Smchd1"),
         col = c("goldenrod","red", "black", "blue", "forestgreen"),
         cex = 1.2, 
         pch = c(15, 16, 16, 16, 17),
         pt.cex = c(1.2, 1, 0.3, 1, 1.5))
dev.off()

# Compare to short-read results from Chen et al.
up <- read.csv("./Chen_et_al_2015_PNAS/Upregulated_in_Smchd1_null.csv", skip=1)
down <- read.csv("./Chen_et_al_2015_PNAS/Downregulated_in_Smchd1_null.csv", skip=1)
index.up <- x$genes$ENTREZID %in% up$EntrezID
index.down <- x$genes$ENTREZID %in% down$EntrezID
index.all <- index.up | index.down
gene.weights <- rep(0, nrow(x))
m <- match(up$EntrezID, x$genes$ENTREZID)
gene.weights[m[!is.na(m)]] <- up$logFC[!is.na(m)]
m <- match(down$EntrezID, x$genes$ENTREZID)
gene.weights[m[!is.na(m)]] <- down$logFC[!is.na(m)]

up.strict <- up[up$adj.P.Val < 0.0001, ]
down.strict <- down[down$adj.P.Val < 0.0001, ]
dim(up.strict)
dim(down.strict)
index.up.strict <- x$genes$ENTREZID %in% up.strict$EntrezID
index.down.strict <- x$genes$ENTREZID %in% down.strict$EntrezID
index.strict <- index.up.strict | index.down.strict
# calculate gene weight
gene.weights <- rep(0, nrow(x))
m <- match(up.strict$EntrezID, x$genes$ENTREZID)
gene.weights[m[!is.na(m)]] <- up.strict$logFC[!is.na(m)]
m <- match(down.strict$EntrezID, x$genes$ENTREZID)
gene.weights[m[!is.na(m)]] <- down.strict$logFC[!is.na(m)]

print(roast(v, index.strict, gene.weights = gene.weights[index.strict], 
            design=design, nrot = 9999))
# Fig 2F
pdf("plots/barcodeStrict.pdf", height = 5, width = 8)
barcodeplot(efit$t[, 2], index=index.strict,
            gene.weights = gene.weights[index.strict],
            labels = c("WT", "Smchd1-null"),
            xlab = "Moderated t",
            cex.lab = 1.8, cex.axis = 2, cex = 1.5)
# legend("top", legend = c("Short-read up genes", "Short-read down genes"),
#        col = c("red", "blue"),
#          text.col = c("red", "blue"),
#          cex = 1.2)s
dev.off()

# Compare t-statistic and logFC
tt <- topTable(efit, number=Inf)
tt.chen <- read.table("./Chen_et_al_2015_PNAS/TopTable.txt", sep="\t")
m <- match(tt$ENTREZID, tt.chen$EntrezID)
cor(tt$logFC, -tt.chen$logFC[m], use = "complete.obs")
cor(tt$t, -tt.chen$t[m], use = "complete.obs")
m.up <- match(up.strict$EntrezID, tt$ENTREZID)
m.down <- match(down.strict$EntrezID, tt$ENTREZID)

# Supp fig S5
library(cowplot)
pdf("plots/LongvsShorttAndLogFC.pdf", height = 5, width = 8)
par(mfrow=c(1,2))
plot.fc <- smoothScatter(tt$logFC, -tt.chen$logFC[m],
              xlab = "long read log FC",
              ylab = "short read log FC",
              )
abline(coef = c(0, 1), lty = 2)
plot.t <- smoothScatter(tt$t, -tt.chen$t[m],
              xlab = "long read t statistic",
              ylab = "short read t statistic",
              )
abline(coef = c(0, 1),  lty=2)
dev.off()

sessionInfo()