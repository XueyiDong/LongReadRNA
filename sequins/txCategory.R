# setwd("~/analysis/2020/smchd1/sequins")
# flames <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/flames/transcript_count.csv", stringsAsFactors = FALSE)
# flair <- read.delim("./flair/counts_matrix.tsv", stringsAsFactors = FALSE)
# talon <- read.delim("./talon/sequins_talon_talon_abundance_filtered.tsv", stringsAsFactors = FALSE)
# anno <- read.delim("~/annotation/sequins/rnasequin_isoforms_2.4.tsv", stringsAsFactors = FALSE)
# 
# in_both <- c(
#   length(intersect(flames$transcript_id, anno$NAME)),
#   length(intersect(flair$ids, anno$NAME)),
#   length(intersect(talon$annot_transcript_id, anno$NAME))
# )
# 
# # number of novel isoforms
# not_in_ref <- c(
#   as.vector(table(is.na(match(flames$transcript_id, anno$NAME))))[2],
#   as.vector(table(is.na(match(flair$ids, anno$NAME))))[1],
#   as.vector(table(is.na(match(talon$annot_transcript_id, anno$NAME))))[2]
# ) 
# 
# # number of missed isoforms
# not_in_results <- c(
#   as.vector(table(is.na(match(anno$NAME, flames$transcript_id))))[2],
#   as.vector(table(is.na(match(anno$NAME, flair$ids))))[1],
#   as.vector(table(is.na(match(anno$NAME, talon$annot_transcript_id))))[2]
# ) 
# 
# stats <- data.frame(
#   method = rep(c("flames", "FLAIR", "TALON"), 3),
#   class = rep(c("in_both", "not_in_ref", "not_in_results"), c(3, 3, 3)),
#   number = c(in_both, not_in_ref, not_in_results)
# )
# stats$class <- factor(stats$class, levels=c("not_in_ref", "not_in_results", "in_both"))
# 
# library(ggplot2)
# 
# pdf("category.pdf", height = 5, width = 8)
# ggplot(stats, aes(x = method, y=number, fill=class)) +
#   geom_bar(stat="identity") +
#   theme_bw() +
#   labs(y = "number of transcript") +
#   scale_fill_manual(values=c("#5B5C63","#3A668C","#C8819D"))
# dev.off()

# plot flames isoform classification by SQANTI
# read SQANTI class
isoClass.flames <- read.delim("./flames/SQANTI/isoform_annotated.filtered_classification.txt", stringsAsFactors = FALSE)
isoClass.flair <- read.delim("./flair/SQANTI/flair.collapse.isoforms_classification.txt", stringsAsFactors = FALSE)
isoClass.talon <- read.delim("./talon/SQANTI/sequins_talon_talon_classification.txt", stringsAsFactors = FALSE)
isoClass.flames$method <- "FLAMES"
isoClass.flair$method <- "FLAIR"
isoClass.talon$method <- "TALON"
# read count matrix
count.flames <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/flames/transcript_count.csv", stringsAsFactors = FALSE)
count.flair <- read.delim("./flair/counts_matrix.tsv", stringsAsFactors = FALSE)
count.talon <- read.delim("./talon/sequins_talon_talon_abundance_filtered.tsv", stringsAsFactors = FALSE)
count.talon_unfiltered <- read.delim("./talon/sequins_talon_talon_abundance.tsv", stringsAsFactors = FALSE)

isoClass.flames$count <- rowSums(count.flames[, 3:6])[match(isoClass.flames$isoform, paste0("transcript:", count.flames$transcript_id))]
isoClass.flair$count <- rowSums(count.flair[, 2:5])[match(substring(isoClass.flair$isoform, 1, 36), substring(count.flair$ids, 1, 36))]
isoClass.talon$count <- rowSums(count.talon[, 12:15])[match(isoClass.talon$isoform, count.talon$annot_transcript_id)]

isoClass <- rbind(isoClass.flames, isoClass.flair, isoClass.talon)
isoClass$structural_category <- factor(isoClass$structural_category, levels =c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "genic_intron", "intergenic", "fusion", "antisense"))
# isoClass$method <- gsub("flames", "FLAMES", isoClass$method)

library(RColorBrewer)
library(scales)
library(ggplot2)
library(cowplot)
col.category <- brewer.pal(nlevels(isoClass$structural_category), "Set1")

# pdf("isoformClass.pdf", height = 5, width = 8)
plot_isoformClass <- ggplot(isoClass, aes(x = method, fill=structural_category)) + 
  geom_bar(position = position_stack(reverse = TRUE)) +
  theme_bw() +
  geom_hline(yintercept = 164, linetype="dashed") +
  labs(y = "Number of transcripts", x = "Method", fill = "Structural category") +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_fill_manual(values = col.category)
# dev.off()

# pdf("isoformClassCount.pdf", height = 5, width = 8)
plot_isoformCount <- ggplot(isoClass, aes(x=method, y=count, fill=structural_category)) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  labs(y = "Read count", x = "Method", fill = "Structural category") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_fill_manual(values = col.category) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))
# dev.off()

leg <- get_legend(plot_isoformCount)

pdf("plots/isoformClassAndCount.pdf", height = 5, width = 17)
plot_grid(plot_isoformClass, plot_isoformCount + theme(legend.position = "none"), leg, rel_heights = c(1, 0.2))
dev.off()

pdf("IsoformClassLengthflames.pdf")
ggplot(isoClass.flames, aes(x=structural_category, y=length, fill=structural_category)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  scale_fill_manual(values = col.category[1:4]) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        text = element_text(size = 16), legend.position = "none") +
  labs(x = "Structural category", y = "Length")
dev.off()

plot(isoClass$diff_to_TSS, isoClass$diff_to_TTS)
