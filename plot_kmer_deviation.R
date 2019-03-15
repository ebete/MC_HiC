#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
  library(ggpubr)
  library(dplyr)
})


reference_2mer <- read.csv("/data0/thom/mm9/mm9_2mer.csv", "\t", sep = "\t", header = T) %>%
  mutate(rel_count = count / sum(count))
window_2mer <- read.csv("/data0/thom/mc4c_window/mc4c_window_2mer.csv", "\t", sep = "\t", header = T) %>%
  mutate(rel_count = count / sum(count))
mapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/mappable_2mer.csv", "\t", sep = "\t", header = T) %>%
  mutate(rel_count = count / sum(count))
unmapped_2mer <- read.csv("/data0/thom/mergemap_splitmap/unmapped_2mer.csv", sep = "\t", header = T) %>%
  mutate(rel_count = count / sum(count))


counts_2mer <- data.frame(
  first = mapped_2mer$first,
  second = mapped_2mer$second,
  seq = as.factor(paste0(as.character(mapped_2mer$first), as.character(mapped_2mer$second))),
  mapped = mapped_2mer$rel_count,
  unmapped = unmapped_2mer$rel_count,
  reference = reference_2mer$rel_count,
  window = window_2mer$rel_count
)
counts_2mer.m <- melt(counts_2mer, id.vars = c("first", "second", "seq", "reference")) %>%
  mutate(deviation_from_expected = value - 4 ^ - 2, deviation_from_reference = value - reference)
# max_count <- max(abs(counts_2mer.m[counts_2mer.m$variable %in% c("unmapped", "mapped"), ]$value))

min_difference <- 0.025
ggplot(counts_2mer.m, aes(x = variable, y = deviation_from_reference, fill = variable)) +
  geom_rect(data = subset(counts_2mer.m, variable == "window"), inherit.aes = F, aes(fill = ifelse(abs(deviation_from_reference) >= min_difference, "2", "1")), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_col(position = "identity", alpha = 0.7) +
  # geom_segment(inherit.aes = F, aes(x = 0, xend = 3, y = reference, yend = reference, color = "reference"), size = 1, linetype = "dotted") +
  facet_grid(second ~ first, labeller = label_both) +
  #  geom_text(data = subset(counts_2mer.m, variable == "reference"), aes(x=1,y=-0.04, label = sprintf("%.1f%%", value*100))) +
  scale_y_continuous(labels = scales::percent, breaks = pretty_breaks()) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_pubr(legend = "top", border = T) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing = unit(0, units = "mm"),
    strip.background = element_rect(fill = "gray98", colour = "black")
  ) +
  scale_fill_brewer(breaks = c("unmapped", "mapped", "window"), palette = "Set1") +
  # scale_color_brewer(breaks = c(), palette = "Dark2") +
  labs(
    x = "",
    y = "Deviation from reference",
    title = "Nucleotide association between mapped and unmapped reads"
  )
