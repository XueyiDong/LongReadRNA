library(ggplot2)

qcdata <- readRDS("./summaryInfo.RDS")
qcdata2 <- readRDS("./summaryInfo2.RDS")
qcdata <- as.data.frame(qcdata, stringsAsFactors = FALSE)
qcdata2 <- as.data.frame(qcdata2, stringsAsFactors = FALSE)
qcdata <- rbind(qcdata, qcdata2)
rm(qcdata2)
qcdata$Read_length <- as.numeric(qcdata$Read_length)
qcdata$Qscore <- as.numeric(qcdata$Qscore)

maxLength = max(qcdata$Read_length)
maxLength
qcdata$LengthGroup <- Hmisc::cut2(qcdata$Read_length, cuts = c(0, 500, 1000,
                                                     3000, maxLength))
qcdata$LengthGroup <- gsub(" ", "", qcdata$LengthGroup)
qcdata$LengthGroup <- gsub(",", ", ", qcdata$LengthGroup)

qcdata$LengthGroup <- factor(qcdata$LengthGroup, levels = c(
  "[0, 500)", "[500, 1000)", "[1000, 3000)", "[3000, 1353320]"
))
# Fig 1C
pdf("plots/lengthQualityViolin.pdf", height = 4, width = 8)
ggplot(qcdata, aes(x=LengthGroup, y=Qscore, fill=LengthGroup, colour = LengthGroup)) +
  geom_violin(alpha = 0.4) +
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none") +
  labs(x = "Read length")
dev.off()

