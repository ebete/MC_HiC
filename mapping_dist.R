#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
})

# draw mapped fragments on top of the original read
aln_coords <- read.csv("/data0/thom/conservative_aln/aln_coords.csv", sep = ";", header = T)
plot.new(); plot.window(ylim = c(1, nrow(aln_coords) + 1), xlim = c(min(aln_coords$start) - 1, max(aln_coords$end) + 1))
for (i in 1 : nrow(aln_coords)) {
  segments(aln_coords[i,]$start, i, aln_coords[i,]$end, i, col = aln_coords[i,]$ref, lty = aln_coords[i,]$strand + 1, lwd = 2)
}

# plot the distribution of the distance to DpnII sites
cutsite_dist <- read.csv("/data0/thom/conservative_aln/LVR_HS5_NP-digested_chimeric.csv", sep = "\t", header = T)
cutsite_dist.melt <- melt(cutsite_dist)
ggplot(cutsite_dist.melt, aes(x = - value, fill = variable)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 10) +
  scale_x_continuous(limits = c(- 200, 200)) +
  scale_y_continuous() +
#  coord_cartesian(ylim = c(0, 1e5)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  xlab("Distance to nearest DpnII from fragment end") +
  ggtitle("Distribution of distance between mapped fragments and DpnII")

# plot the distribution of mapped fragment lengths
map_len <- data.frame(x = scan("/data0/thom/conservative_aln/LVR_HS5_NP-digested_chimeric_len.txt", numeric()))
ggplot(map_len, aes(x = x, fill = "Mapped length")) +
  geom_density() +
  scale_x_continuous(limits = c(0, 200)) +
  scale_y_continuous() +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()
