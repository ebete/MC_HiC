#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  library(ggbio)
  library(gridExtra)
  library(cowplot)
})

interactions <- read.csv("/data0/thom/temp/bwa_fixed.bam.csv_interactions.csv", sep = ";")
chr7 <- interactions[(interactions$chr.1 == "chr7") & (interactions$chr.2 == "chr7"),]

# interaction map
g <- ggplot(chr7, aes(x = pos.1, y = pos.2)) +
  theme_minimal() +
  scale_y_reverse(labels = function(x) { sprintf("%.0f Mb", x / 1e6)}, limits = c(115e6, 105e6)) +
  scale_x_continuous(position = "top", labels = function(x) { sprintf("%.0f Mb", x / 1e6)}, limits = c(105e6, 115e6)) +
  xlab("chr7") +
  ylab("chr7") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
plt1 <- g +
  geom_bin2d(bins = 500) +
  scale_fill_continuous(type = "viridis", limits = c(0, 30))

# reference genome annotation
genemap <- autoplot(TxDb.Mmusculus.UCSC.mm9.knownGene, which = range(GRanges(Rle(c("chr7"), c(1)), IRanges(105e6, width = 10e6))), gap.geom = "chevron", label = F)
plt2 <- genemap@ggplot +
  ggtitle("chr7 MC-4C interaction map (105 Mb - 115 Mb)") +
  scale_x_continuous(labels = function(x) { sprintf("%.0f Mb", x / 1e6)}) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# plot both
plot_grid(plt2, plt1, ncol = 1, align = "h", axis = "t", rel_heights = c(1 / 5, 4 / 5))
