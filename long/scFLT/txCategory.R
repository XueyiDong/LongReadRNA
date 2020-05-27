setwd("./analysis/smchd1/long/scFLT")
library(ggplot2)

isoClass <- read.delim("isoform_annotated_classification.txt", stringsAsFactors = FALSE)
isoClass$structural_category <- factor(isoClass$structural_category, levels =c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "genic_intron", "fusion", "antisense"))

pdf("txCategory.pdf", height = 5, width = 8)
ggplot(isoClass, aes(x=structural_category, fill=structural_category)) +
  geom_bar() +
  labs(y = "Number of transcripts") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

counts <- read.csv("./run1/transcript_count.csv")
isoClass$count <- rowSums(counts[3:16])[match(isoClass$isoform, paste0("transcript:", counts$transcript_id))]

pdf("txCategoryCount.pdf", height = 5, width = 8)
ggplot(isoClass[!is.na(isoClass$count),], aes(x=structural_category, y=count, fill=structural_category))+
  geom_bar(stat="identity")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

classCount <- aggregate(isoClass$count[!is.na(isoClass$count)], by=list(structural_category=isoClass$structural_category[!is.na(isoClass$count)]), FUN=sum)
