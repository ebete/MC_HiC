#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
})

unmapped <- scan("/data0/thom/pipeline/unmapped_len.txt", numeric())
mapped <- scan("/data0/thom/pipeline/mapped_len.txt", numeric())

df <- rbind(
data.frame(source = "cis mapped", length = mapped),
data.frame(source = "trans mapped", length = unmapped)
)

ggplot(df, aes(x = length, fill = source)) +
  geom_density(alpha = 0.5) +
#  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous() +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold"),
  legend.position = "bottom"
  ) +
  xlab("Read length") +
  ylab("Density") +
  ggtitle("Distribution of read lengths")
