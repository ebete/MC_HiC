#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
})


mapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/mappable_2mer.csv", "\t", sep = "\t", header = T)
unmapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/unmapped_2mer.csv", sep = "\t", header = T)

mapped_2mer$count <- mapped_2mer$count / sum(mapped_2mer$count)
unmapped_2mer$count <- unmapped_2mer$count / sum(unmapped_2mer$count)

counts_2mer <- data.frame(seq = mapped_2mer$X2mer, mapped = mapped_2mer$count, unmapped = unmapped_2mer$count)
counts_2mer.m <- melt(counts_2mer, id = "seq")

counts_2mer <- data.frame(seq = mapped_2mer$X2mer, mapped = mapped_2mer$count, unmapped = unmapped_2mer$count)
counts_2mer.m <- melt(counts_2mer, id = "seq")
counts_2mer.m$value <- counts_2mer.m$value - 1 / 4 ^ 2

ggplot(counts_2mer.m, aes(x = seq, y = value, fill = variable)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(- max(abs(counts_2mer.m$value)), max(abs(counts_2mer.m$value)))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
  title = "2-mer rate difference betweeen mappable and unmappable reads",
  x = "2-mer",
  y = "Deviation from expected (1 / 4^2)"
  )
