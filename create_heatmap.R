#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
  library(gplots)
  library(RColorBrewer)
})

# functions
relfreq <- function(x) { x / sum(x)}

# load data and normalise
df <- read.csv("/home/thom/mc_hic/mc_4c/demux/result.pdf.csv")[, - 1]
mtx <- t(data.matrix(df))
mtx.norm <- t(apply(mtx, c(1), relfreq))

# plot dendrogram
col_palette <- colorRampPalette(brewer.pal(n = 9, name = "Reds"))
heatmap.2(mtx.norm,
dendrogram = "row",
Colv = F,
margins = c(4, 12),
col = col_palette(100),
tracecol = "gray",
cellnote = mtx,
notecol = "black"
)
