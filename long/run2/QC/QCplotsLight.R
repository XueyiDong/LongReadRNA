library(ggplot2)
readnum <- data.frame(
  number = c(19489349, 17668058, 13854515, 22164601, 6825392, 2657225, 18833786, 17433027, 13854516, 20677451, 5986242, 2686678),
  run = c(rep("run 1", 6), rep("run 2", 6)),
  quality = rep(rep(c("Pass", "Fail"), c(3, 3)), 2),
  category = factor(c(rep(c("raw", "demultiplexed",  "assigned"), 4)), levels= c("raw", "demultiplexed",  "assigned"))
)
pdf("totalreadnum.pdf", height = 4, width = 8)
ggplot(readnum, aes(x=category, y=number, fill=quality)) +
  geom_bar(stat = "identity") +
  facet_grid(.~run) +
  theme_bw() +
  geom_text(aes(label=number), position = position_stack(vjust = 0.5), size = 3)
dev.off()

readnum.ag <- aggregate(x=readnum$number, by=list(readnum$quality, readnum$category), FUN=sum)
colnames(readnum.ag) <- c("quality", "category", "number")
pdf("totalreadnum_combinedBatch.pdf", height = 4, width = 8)
ggplot(readnum.ag, aes(x=category, y=number, fill=quality)) +
  geom_bar(stat = "identity") +
  # facet_grid(.~run) +
  theme_bw() +
  geom_text(aes(label=number), position = position_stack(vjust = 0.5), size = 3)
dev.off()
