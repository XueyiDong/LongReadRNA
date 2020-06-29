setwd("~/analysis/smchd1/long/sequins")
fltsa <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/transcript_count.csv", stringsAsFactors = FALSE)
flair <- read.delim("./flair/counts_matrix.tsv", stringsAsFactors = FALSE)
talon <- read.delim("./talon/sequins_talon_talon_abundance_filtered.tsv", stringsAsFactors = FALSE)
anno <- read.delim("~/annotation/sequins/rnasequin_isoforms_2.4.tsv", stringsAsFactors = FALSE)

in_both <- c(
  length(intersect(fltsa$transcript_id, anno$NAME)),
  length(intersect(flair$ids, anno$NAME)),
  length(intersect(talon$annot_transcript_id, anno$NAME))
)

# number of novel isoforms
not_in_ref <- c(
  as.vector(table(is.na(match(fltsa$transcript_id, anno$NAME))))[2],
  as.vector(table(is.na(match(flair$ids, anno$NAME))))[1],
  as.vector(table(is.na(match(talon$annot_transcript_id, anno$NAME))))[2]
) 

# number of missed isoforms
not_in_results <- c(
  as.vector(table(is.na(match(anno$NAME, fltsa$transcript_id))))[2],
  as.vector(table(is.na(match(anno$NAME, flair$ids))))[1],
  as.vector(table(is.na(match(anno$NAME, talon$annot_transcript_id))))[2]
) 

stats <- data.frame(
  method = rep(c("FLTSA", "FLAIR", "TALON"), 3),
  class = rep(c("in_both", "not_in_ref", "not_in_results"), c(3, 3, 3)),
  number = c(in_both, not_in_ref, not_in_results)
)
stats$class <- factor(stats$class, levels=c("not_in_ref", "not_in_results", "in_both"))

library(ggplot2)

pdf("category.pdf", height = 5, width = 8)
ggplot(stats, aes(x = method, y=number, fill=class)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(y = "number of transcript") +
  scale_fill_manual(values=c("#5B5C63","#3A668C","#C8819D"))
dev.off()

# plot FLTSA isoform classification by SQANTI
# read SAANTI class
isoClass.fltsa <- read.delim("./SQANTI/isoform_annotated_classification.txt", stringsAsFactors = FALSE)
isoClass.flair <- read.delim("./flair/SQANTI/flair.collapse.rmScBc.isoforms_classification.txt", stringsAsFactors = FALSE)
isoClass.talon <- read.delim("./talon/SQANTI/sequins_talon_talon_classification.txt", stringsAsFactors = FALSE)
isoClass.fltsa$method <- "FLTSA"
isoClass.flair$method <- "FLAIR"
isoClass.talon$method <- "TALON"
# read count matrix
count.fltsa <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/transcript_count.csv", stringsAsFactors = FALSE)
count.talon <- read.delim("./talon/sequins_talon_talon_abundance_filtered.tsv", sep = "\t", stringsAsFactors = FALSE)
count.talon_unfiltered <- read.delim("./talon/sequins_talon_talon_abundance.tsv", sep = "\t", stringsAsFactors = FALSE)
count.flair <- read.delim("./flair/counts_matrix.tsv", sep = "\t", stringsAsFactors = FALSE)

isoClass.fltsa$count <- rowSums(count.fltsa[, 3:6])[match(isoClass.fltsa$isoform, paste0("transcript:", count.fltsa$transcript_id))]
isoClass.flair$count <- rowSums(count.flair[, 2:5])[match(substring(isoClass.flair$isoform, 1, 36), substring(count.flair$ids, 1, 36))]
isoClass.talon$count <- rowSums(count.talon[, 12:15])[match(isoClass.talon$isoform, count.talon$annot_transcript_id)]

isoClass <- rbind(isoClass.fltsa, isoClass.flair, isoClass.talon)
isoClass$structural_category <- factor(isoClass$structural_category, levels =c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "genic_intron", "intergenic", "fusion", "antisense"))
isoClass$method <- gsub("FLTSA", "FLAMES", isoClass$method)

library(RColorBrewer)
library(scales)
col.category <- brewer.pal(nlevels(isoClass$structural_category), "Set1")

pdf("isoformClass.pdf", height = 5, width = 8)
ggplot(isoClass, aes(x = method, fill=structural_category)) + 
  geom_bar(position = position_stack(reverse = TRUE)) +
  theme_bw() +
  geom_hline(yintercept = 164, linetype="dashed") +
  labs(y = "Number of transcripts", x = "Method", fill = "Structural category") +
  theme(text = element_text(size = 20), legend.position = c(0.75, 0.75),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.background = element_rect(linetype = "solid", colour = "black", size = 0.25)) +
  scale_fill_manual(values = col.category[1:6])
dev.off()

pdf("isoformClassCount.pdf", height = 5, width = 8)
ggplot(isoClass, aes(x=method, y=count, fill=structural_category)) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  labs(y = "Read count", x = "Method", fill = "Structural category") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_fill_manual(values = col.category[1:6]) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))
dev.off()

pdf("IsoformClassLengthFLTSA.pdf")
ggplot(isoClass.fltsa, aes(x=structural_category, y=length, fill=structural_category)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  scale_fill_manual(values = col.category[1:4]) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        text = element_text(size = 16), legend.position = "none") +
  labs(x = "Structural category", y = "Length")
dev.off()

plot(isoClass$diff_to_TSS, isoClass$diff_to_TTS)
