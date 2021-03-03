library(ggplot2)
library(scales)

readnum <- data.frame(
  number = c(81165187, 57234663, 51135000, 37897462, 33311549),
  category = factor(c("raw reads", "quality filtered reads", "demultiplexed reads", "reads from chosen samples", "gene counts"), 
                    levels = c("raw reads", "quality filtered reads", "demultiplexed reads", "reads from chosen samples", "gene counts"))
)
# Fig 1B
pdf("totalreadnum.pdf", height = 4, width = 8)
ggplot(readnum, aes(x=category, y=number)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),) +
  labs(x = "Category", y = "Number") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))
dev.off()
