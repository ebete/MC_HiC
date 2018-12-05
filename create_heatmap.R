#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
  library(gplots)
  library(RColorBrewer)
})

# export to PDF
pdf("create_heatmap.pdf", width = 11.69, height = 8.27, pointsize = 6)

# functions
relfreq <- function(x) { x / sum(x)}

# load data and normalise
df <- read.csv("/home/thom/mc_hic/mc_4c_digested/demux/result.pdf.csv")[- c(1, 2), - 1]
mtx <- t(data.matrix(df))
#mtx.norm <- t(apply(mtx, c(1), relfreq))

# plot dendrogram
col_palette <- colorRampPalette(brewer.pal(n = 9, name = "Reds"))
col_palette <- colorRampPalette(c("red", "white", "green"))

heatmap.2(mtx,
dendrogram = "row",
hclustfun = function(x) {hclust(x, method = "average")},
Colv = F,
margins = c(4, 12),
col = col_palette(15),
tracecol = NULL,
cellnote = mtx,
notecol = "black"
)

dev.off()
