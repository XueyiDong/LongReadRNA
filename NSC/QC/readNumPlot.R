library(ggplot2)
library(scales)
# readnum <- data.frame(
#   number = c(19489349, 17668058, 13854515, 22164601, 6825392, 2657225, 18833786, 17433027, 13854516, 20677451, 5986242, 2686678),
#   run = c(rep("run 1", 6), rep("run 2", 6)),
#   quality = rep(rep(c("Pass", "Fail"), c(3, 3)), 2),
#   category = factor(c(rep(c("raw", "demultiplexed",  "assigned"), 4)), levels= c("raw", "demultiplexed",  "assigned"))
# )
# pdf("totalreadnum.pdf", height = 4, width = 8)
# ggplot(readnum, aes(x=category, y=number, fill=quality)) +
#   geom_bar(stat = "identity") +
#   facet_grid(.~run) +
#   theme_bw() +
#   geom_text(aes(label=number), position = position_stack(vjust = 0.5), size = 3)
# dev.off()
# 
# readnum.ag <- aggregate(x=readnum$number, by=list(readnum$quality, readnum$category), FUN=sum)
# colnames(readnum.ag) <- c("quality", "category", "number")
# pdf("totalreadnum_combinedBatch.pdf", height = 4, width = 8)
# ggplot(readnum.ag, aes(x=category, y=number, fill=quality)) +
#   geom_bar(stat = "identity") +
#   # facet_grid(.~run) +
#   theme_bw() +
#   geom_text(aes(label=number), position = position_stack(vjust = 0.5), size = 5) +
#   theme(text = element_text(size=20)) +
#   labs(x = "", y = "Number of reads")
# dev.off()


readnum <- data.frame(
  number = c(81165187, 57234663, 51135000, 37897462, 33311549),
  category = factor(c("raw reads", "quality filtered reads", "demultiplexed reads", "reads from chosen samples", "gene counts"), 
                    levels = c("raw reads", "quality filtered reads", "demultiplexed reads", "reads from chosen samples", "gene counts"))
)

pdf("totalreadnum.pdf", height = 4, width = 8)
ggplot(readnum, aes(x=category, y=number)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),) +
  labs(x = "Category", y = "Number") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))
dev.off()
