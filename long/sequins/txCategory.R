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
ggplot(stats, aes(x = method, y=number, fill=class)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(y = "number of transcript") +
  scale_fill_manual(values=c("#2A557B","#BB6C8B","#474A52"))
