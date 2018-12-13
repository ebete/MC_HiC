#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  library(ggbio)
  library(gridExtra)
  library(cowplot)
})

interactions <- read.csv("~/mc_hic/mc_4c/bwa_gapped_14seed.bam_interactions.csv", sep = ";")

# interaction map
g <- ggplot(interactions, aes(x = x, y = y)) +
  theme_minimal() +
  scale_y_reverse(labels = scales::scientific) +
  scale_x_continuous(position = "top", labels = scales::scientific) +
  xlab("chr7") +
  ylab("chr7") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
plt1 <- g +
  geom_bin2d(bins = 500) +
  scale_fill_continuous(type = "viridis", limits = c(0, 30))

# reference genome annotation
genemap <- autoplot(TxDb.Mmusculus.UCSC.mm9.knownGene, which = range(GRanges(Rle(c("chr7"), c(1)), IRanges(105000000, width = 10000000))), gap.geom = "chevron", label = F)
plt2 <- genemap@ggplot +
  ggtitle("chr7 MC-4C interaction map (105 Mb - 115 Mb)") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# plot both
plot_grid(plt2, plt1, ncol = 1, align = "h", axis = "t", rel_heights = c(1 / 5, 4 / 5))
