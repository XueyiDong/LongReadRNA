tracking <- read.table("scFLT.tracking", stringsAsFactors = FALSE)
library(ggplot2)
pdf("compare_category.pdf", height = 4)
ggplot(tracking, aes(x = V4, fill=V4)) +
  geom_bar() +
  geom_text(stat="count", aes(label=..count..), vjust=-1) +
  theme_bw()
dev.off()
