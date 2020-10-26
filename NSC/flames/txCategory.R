# setwd("./analysis/smchd1/long/scFLT")
library(ggplot2)
library(RColorBrewer)
library(scales)

isoClass <- read.delim("isoform_annotated.filtered_classification.txt", stringsAsFactors = FALSE)
isoClass$structural_category <- factor(isoClass$structural_category, levels =c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "genic_intron", "intergenic", "fusion", "antisense"))

col.category <- brewer.pal(nlevels(isoClass$structural_category), "Set1")

pdf("txCategory.pdf", height = 5, width = 8)
ggplot(isoClass, aes(x=structural_category, fill=structural_category)) +
  geom_bar() +
  labs(x = "Structural category", y = "Number of transcripts") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), 
        text = element_text(size = 20),
        legend.position = "none") +
  scale_fill_manual(values = col.category[c(1:4, 7:9, 6)])
dev.off()

counts <- read.csv("./results/transcript_count.csv")
isoClass$count <- rowSums(counts[3:10])[match(isoClass$isoform, paste0("transcript:", counts$transcript_id))]

pdf("txCategoryCount.pdf", height = 5, width = 8)
ggplot(isoClass[!is.na(isoClass$count),], aes(x=structural_category, y=count, fill=structural_category))+
  geom_bar(stat="identity")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), 
        text = element_text(size = 20),
        legend.position = "none")+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
  scale_fill_manual(values = col.category[c(1:4, 7:9, 6)]) +
  labs(x = "Structural category", y = "Read count")
dev.off()

classCount <- aggregate(isoClass$count[!is.na(isoClass$count)], by=list(structural_category=isoClass$structural_category[!is.na(isoClass$count)]), FUN=sum)
