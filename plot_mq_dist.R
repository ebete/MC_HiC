#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
})

# load MAPQ matrix
mq <- read.csv("/home/thom/mc_hic/mc_4c/mq_dist.txt", header = F, sep = ";", as.is = T)
mq <- t(mq)
colnames(mq) <- mq[1,]
mq <- mq[- 1,]
class(mq) <- "numeric"
mq.melt <- melt(mq)

# plot density
ggplot(mq.melt, aes(x = value, fill = Var2)) +
  geom_histogram(position = "dodge", na.rm = T, binwidth = 10) +
  theme_bw() +
#  scale_x_continuous(limits = c(0, 100), breaks = 1 : 10 * 10) +
  scale_y_continuous(labels = comma) +
  xlab("MAPQ") +
  ylab("# reads") +
  ggtitle("Distribution of MAPQ values") +
  theme(
  plot.title = element_text(face = "bold", hjust = 0.5)
  )
