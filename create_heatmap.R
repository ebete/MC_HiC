#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
  library(gplots)
  library(RColorBrewer)
  library(abind)
})

#par(bg = 'black', fg = 'white')
# functions
relfreq <- function(x) { x / sum(x)}

# load data
df <- read.csv("/home/thom/mc_hic/mc_4c_digested/demux/result.pdf.csv")[- 1 : - 3, - 1]
mtx <- t(data.matrix(df))

# infer score
mtx.scores <- matrix(mtx[, 1])
row.names(mtx.scores) <- row.names(mtx)
for (i in 2 : ncol(mtx)) {
  mtx.scores <- mtx.scores + mtx[, i] * i
}
# sort based on score
mtx <- mtx[order(mtx.scores[, 1], decreasing = T),]

#mtx.norm <- t(apply(mtx, c(1), relfreq))
# color labels based on sample source
mtx.fac <- rep("black", nrow(mtx))
for (i in 1 : nrow(mtx)) {
  source <- strsplit(row.names(mtx)[i], "_", fixed = T)[[1]]
  source_dataset <- source[length(source)]  # get data set
  source_conf <- paste(source[- length(source)], collapse = "_")  # get alignment conf
  mtx.fac[i] <- source_dataset
}
mtx.fac <- as.factor(mtx.fac)
mtx.col <- colorRampPalette(brewer.pal(n = length(levels(mtx.fac)), name = "Set1"))(length(levels(mtx.fac)))[as.numeric(mtx.fac)]

rankings <- data.frame(rank = 1 : 13)
run <- 1
for (l in levels(mtx.fac)) {
  ranked <- data.frame(row.names(mtx[mtx.fac == l,]))
  colnames(ranked) <- sprintf("run_%d", as.integer(run))
  rankings <- cbind(rankings, ranked)
  run <- run + 1
}


# plot dendrogram
col_palette <- colorRampPalette(brewer.pal(n = 9, name = "PuOr"))
#col_palette <- colorRampPalette(c("orange", "white", "cyan"))

pdf("heatmaps.pdf", width = 11.7, height = 8.27, pointsize = 12, colormodel = "cmyk")
for (l in levels(mtx.fac)) {
  heatmap.2(mtx[mtx.fac == l,],
  dendrogram = "none",
  hclustfun = function(x) {hclust(x, method = "average")},
  Colv = F,
  Rowv = F,
  margins = c(4, 20),
  col = col_palette(15),
  tracecol = NULL,
  cellnote = mtx[mtx.fac == l,],
  notecol = "black",
  #    colRow = mtx.col[mtx.fac == l],
  colCol = "black"
  )
}
dev.off()
