library(ggplot2)
readnum <- data.frame(
  number = c(41653950, 28539637, 25074894, 16511740, 39511237, 27493223, 24919069, 16541194),
  run = c(rep("1", 4), rep("2", 4)),
  category = factor(c(rep(c("raw", "demultiplexed", "mapped", "assigned"))), levels= c("raw", "demultiplexed", "mapped", "assigned"))
)
pdf("totalreadnum.pdf", height = 5)
ggplot(readnum, aes(x=category, y=number, fill=run)) +
  geom_bar(stat = "identity")
dev.off()
