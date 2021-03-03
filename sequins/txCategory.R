# plot flames isoform classification by SQANTI, fig 3E and 3F
# read SQANTI class
isoClass.flames <- read.delim("./flames/SQANTI/isoform_annotated.filtered_classification.txt", stringsAsFactors = FALSE)
isoClass.flair <- read.delim("./flair/SQANTI/flair.collapse.isoforms_classification.txt", stringsAsFactors = FALSE)
isoClass.talon <- read.delim("./talon/SQANTI/sequins_talon_talon_classification.txt", stringsAsFactors = FALSE)
isoClass.flames$method <- "FLAMES"
isoClass.flair$method <- "FLAIR"
isoClass.talon$method <- "TALON"
# read count matrix
count.flames <- read.csv("./flames/results/transcript_count.csv", stringsAsFactors = FALSE)
count.flair <- read.delim("./flair/counts_matrix.tsv", stringsAsFactors = FALSE)
count.talon <- read.delim("./talon/sequins_talon_talon_abundance_filtered.tsv", stringsAsFactors = FALSE)
count.talon_unfiltered <- read.delim("./talon/sequins_talon_talon_abundance.tsv", stringsAsFactors = FALSE)

isoClass.flames$count <- rowSums(count.flames[, 3:6])[match(isoClass.flames$isoform, paste0("transcript:", count.flames$transcript_id))]
isoClass.flair$count <- rowSums(count.flair[, 2:5])[match(substring(isoClass.flair$isoform, 1, 36), substring(count.flair$ids, 1, 36))]
isoClass.talon$count <- rowSums(count.talon[, 12:15])[match(isoClass.talon$isoform, count.talon$annot_transcript_id)]

isoClass <- rbind(isoClass.flames, isoClass.flair, isoClass.talon)
isoClass$structural_category[isoClass$structural_category %in% c("genic", "genic_intron", "intergenic")] <- "genic/intergenic"
isoClass$structural_category <- factor(isoClass$structural_category, levels =c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic/intergenic", "fusion", "antisense"))

library(RColorBrewer)
library(scales)
library(ggplot2)
library(cowplot)
col.category <- brewer.pal(nlevels(isoClass$structural_category), "Set1")

plot_isoformClass <- ggplot(isoClass, aes(x = method, fill=structural_category)) + 
  geom_bar(position = position_stack(reverse = TRUE)) +
  theme_bw() +
  geom_hline(yintercept = 164, linetype="dashed") +
  labs(y = "Number of transcripts", x = "Method", fill = "Structural category") +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_fill_manual(values = col.category)

plot_isoformCount <- ggplot(isoClass, aes(x=method, y=count, fill=structural_category)) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  labs(y = "Read count", x = "Method", fill = "Structural category") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_fill_manual(values = col.category) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))

leg <- get_legend(plot_isoformCount)

pdf("plots/isoformClassAndCount.pdf", height = 5, width = 17)
plot_grid(plot_isoformClass, plot_isoformCount + theme(legend.position = "none"), leg, rel_heights = c(1, 0.2))
dev.off()