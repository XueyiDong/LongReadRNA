## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------
library(ggplot2)


## --------------------------------------------------------------------------
qcdata <- readRDS("./summaryInfo.RDS")


## --------------------------------------------------------------------------
pdf("plots/length.pdf", height = 5)
ggplot(qcdata, aes(x = Read_length )) +
  geom_density() +
  theme_bw() +
  coord_trans(x = "log10")
dev.off()

pdf("plots/lengthSample.pdf", height = 5)
ggplot(qcdata, aes(x = Read_length, colour=Barcode )) +
  geom_density() +
  theme_bw() +
  coord_trans(x = "log10")
dev.off()


## --------------------------------------------------------------------------
cutoff <- quantile(qcdata$Read_length, 0.99)
cat("cutoff (read length 0.99 quantile): ", cutoff, "\n")



## --------------------------------------------------------------------------
# pdf("lengthQscore.pdf", )
# ggplot(qcdata, aes(x = Read_length, y = Qscore)) +
#   geom_point() +
#   xlim(NA, cutoff)
# dev.off()


## --------------------------------------------------------------------------
pdf("plots/lengthpassfail.pdf", height = 4, width = 8)
qcdata$pass <- qcdata$Qscore > 7
ggplot(qcdata, aes(x = Read_length, fill = pass)) +
  geom_density(alpha = 0.4) +
  xlim(NA, cutoff) +
  theme_bw()
dev.off()


## --------------------------------------------------------------------------
cat("1/3: ", quantile(qcdata$Read_length, 1/3))
cat("2/3: ", quantile(qcdata$Read_length, 2/3))
cat("0.95: ", quantile(qcdata$Read_length, 0.95))


## --------------------------------------------------------------------------
qcdata$length <- "medium"
qcdata$length[qcdata$Read_length < quantile(qcdata$Read_length, 1/3)] <- "short"
qcdata$length[qcdata$Read_length > quantile(qcdata$Read_length, 2/3)] <- "long"
qcdata$length[qcdata$Read_length > quantile(qcdata$Read_length, 0.95)] <- "extra long"
pdf("plots/lengthQualityviolin.pdf", height = 4, width = 8)
ggplot(qcdata, aes(x = length, y = Qscore, fill = length, colour = length)) +
  geom_violin(alpha = 0.4) +
  theme_bw() +
  theme(text = element_text(size=20))
dev.off()


## --------------------------------------------------------------------------
maxLength = max(qcdata$Read_length)
qcdata$LengthGroup <- Hmisc::cut2(qcdata$Read_length, cuts = c(0, 500, 1000,
                                                     3000, maxLength))
qcdata$LengthGroup <- gsub(" ", "", qcdata$LengthGroup)
qcdata$LengthGroup <- gsub(",", ", ", qcdata$LengthGroup)

qcdata$LengthGroup <- factor(qcdata$LengthGroup, levels = c(
  "[0, 500)", "[500, 1000)", "[1000, 3000)", "[3000, 1823562]"
))

pdf("plots/lengthGroupQualityViolin.pdf", height = 4, width = 8)
ggplot(qcdata, aes(x=LengthGroup, y=Qscore, fill=LengthGroup, colour = LengthGroup)) +
  geom_violin(alpha = 0.4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none") +
  labs(x = "Read length")
dev.off()

