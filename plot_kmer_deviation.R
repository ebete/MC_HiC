#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
  library(ggpubr)
  library(dplyr)
})


reference_2mer <- read.csv("/data0/thom/mc4c_window/mc4c_window_2mer.csv", "\t", sep = "\t", header = T) %>%
mutate(rel_count = count / sum(count))
mapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/mappable_2mer.csv", "\t", sep = "\t", header = T) %>%
mutate(rel_count = count / sum(count) - reference_2mer$rel_count)
unmapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/unmapped_2mer.csv", sep = "\t", header = T) %>%
mutate(rel_count = count / sum(count) - reference_2mer$rel_count)

counts_2mer <- data.frame(
first = mapped_2mer$first,
second = mapped_2mer$second,
seq = as.factor(paste0(as.character(mapped_2mer$first), as.character(mapped_2mer$second))),
mapped = mapped_2mer$rel_count,
unmapped = unmapped_2mer$rel_count,
reference = reference_2mer$rel_count
)
counts_2mer.m <- melt(counts_2mer, id.vars = c("first", "second", "seq")) %>%
mutate(deviation_from_expected = value - 4 ^ - 2)
max_count <- max(abs(counts_2mer.m[counts_2mer.m$variable != "reference",]$value))

ggplot(subset(counts_2mer.m, variable != "reference"), aes(x = variable, y = value, fill = variable)) +
  geom_col(position = "identity", alpha = 0.5) +
  facet_grid(second ~ first, labeller = label_both) +
#  geom_text(data = subset(counts_2mer.m, variable == "reference"), aes(x=1,y=-0.04, label = sprintf("%.1f%%", value*100))) +
  scale_y_continuous(labels = scales::percent, breaks = pretty_breaks()) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_pubr(legend = "top") +
  theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank(),
  panel.spacing = unit(0, units = "mm"),
  strip.background = element_rect(fill = "gray98", colour = "black"),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  scale_fill_brewer(palette = "Set1") +
  labs(
  x = "",
  y = "Deviation from reference",
  title = "Nucleotide association between mapped and unmapped reads"
  )
