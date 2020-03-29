tracking <- read.table("comp.tracking", stringsAsFactors = FALSE)
library(ggplot2)
pdf("compare_category_sequins.pdf", height = 4)
ggplot(tracking, aes(x = V4, fill=V4)) +
  geom_bar() +
  geom_text(stat="count", aes(label=..count..), vjust=-0.2) +
  theme_bw()
dev.off()

# barplot by  read number but not number of transcripts

count <-read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/FLTSA_output/transcript_count.csv")

info <- strsplit2(tracking$V5, "|", fixed=TRUE)
tracking$transcript <- sub("transcript:", "", info[,2])
count$all <- rowSums(count[,3:ncol(count)])
count$category <- tracking$V4[match(count$transcript_id, tracking$transcript)]
head(count)
# known <- grepl("^ENSMUST", count$transcript_id)
# count$category[known] <- "="
count[(is.na(count$category)),]
count[(is.na(count$category)), "category"] = "="
table(count$category)
count_cat <- aggregate(count$all, by = list(count$category), FUN=sum)
colnames(count_cat) <- c("category", "count")
pdf("compare_category_count.pdf", height = 4)
ggplot(count_cat, aes(x = category, y=count, fill=category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=count), vjust = -0.1, size = 3) +
  theme_bw()
dev.off()