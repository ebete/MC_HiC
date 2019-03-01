#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
  library(latex2exp)
  library(ggpubr)
})


mapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/mappable_2mer.csv", "\t", sep = "\t", header = T)
unmapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/unmapped_2mer.csv", sep = "\t", header = T)

mapped_2mer$rel_count <- mapped_2mer$count / sum(mapped_2mer$count)
unmapped_2mer$rel_count <- unmapped_2mer$count / sum(unmapped_2mer$count)

counts_2mer <- data.frame(
first = mapped_2mer$first,
second = mapped_2mer$second,
seq = as.factor(paste0(as.character(mapped_2mer$first), as.character(mapped_2mer$second))),
mapped = mapped_2mer$rel_count,
unmapped = unmapped_2mer$rel_count
)
counts_2mer.m <- melt(counts_2mer, id.vars = c("first", "second", "seq"))
counts_2mer.m$deviation <- counts_2mer.m$value - 1 / 4 ^ 2

ggplot(counts_2mer.m, aes(x = seq, y = deviation, fill = variable)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(- max(abs(counts_2mer.m$deviation)), max(abs(counts_2mer.m$deviation)))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
  title = "2-mer rate difference betweeen mappable and unmappable reads",
  x = "2-mer",
  y = TeX("Deviation from expected ($4^{-2}$)")
  )

first <- ggplot(counts_2mer.m, aes(x = second, y = deviation, fill = variable)) +
  geom_col(position = "identity", alpha = 0.5) +
  facet_grid(. ~ first, labeller = label_both) +
  scale_y_continuous(limits = c(- max(abs(counts_2mer.m$deviation)), max(abs(counts_2mer.m$deviation))), labels = scales::percent) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(
  axis.text.x = element_text(angle = 90),
  panel.spacing.x = unit(0, units = "mm"),
  strip.background = element_rect(fill = "gray98", colour = "black"),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  scale_fill_brewer(palette = "Set1") +
  labs(
  x = "Second nucleotide",
  y = TeX("Deviation from expected ($4^{-2}$)")
  )

second <- ggplot(counts_2mer.m, aes(x = first, y = deviation, fill = variable)) +
  geom_col(position = "identity", alpha = 0.5) +
  facet_grid(. ~ second, labeller = label_both) +
  scale_y_continuous(limits = c(- max(abs(counts_2mer.m$deviation)), max(abs(counts_2mer.m$deviation))), labels = scales::percent) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(
  axis.text.x = element_text(angle = 90),
  panel.spacing.x = unit(0, units = "mm"),
  strip.background = element_rect(fill = "gray98", colour = "black"),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  scale_fill_brewer(palette = "Set1") +
  labs(
  x = "First nucleotide",
  y = TeX("Deviation from expected ($4^{-2}$)")
  )

combined <- ggarrange(first, second, nrow = 2, labels = c("(1)", "(2)"), hjust = 0, align = "h", common.legend = T, legend = "top")
annotate_figure(combined, top = "Nucleotide association difference betweeen mappable and unmappable reads")
