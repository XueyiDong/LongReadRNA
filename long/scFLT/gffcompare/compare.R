tracking <- read.table("scFLT.tracking", stringsAsFactors = FALSE)
library(ggplot2)
pdf("compare_category.pdf", height = 4)
ggplot(tracking, aes(x = V4, fill=V4)) +
  geom_bar() +
  geom_text(stat="count", aes(label=..count..), vjust=-0.2) +
  theme_bw()
dev.off()

# barplot by  read number but not number of transcripts

count_file <- "../run1/transcript_count.csv"
count <- read.csv(count_file, stringsAsFactors = FALSE)

info <- strsplit2(tracking$V5, "|", fixed=TRUE)
tracking$transcript <- sub("transcript:", "", info[,2])
count$all <- rowSums(count[,3:16])
count$category <- tracking$V4[match(count$transcript_id, tracking$transcript)]
head(count)
known <- grepl("^ENSMUST", count$transcript_id)
count$category[known] <- "="
table(count$category)
count_cat <- aggregate(count$all, by = list(count$category), FUN=sum)
colnames(count_cat) <- c("category", "count")
pdf("compare_category_count.pdf", height = 4)
ggplot(count_cat, aes(x = category, y=count, fill=category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=count), vjust = -0.1, size = 3) +
  theme_bw()
dev.off()