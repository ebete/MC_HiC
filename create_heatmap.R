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
df <- read.csv("/home/thom/mc_hic/mc_4c_digested/demux/result.pdf.csv")[3 : 8, - 1]
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
  source <- strsplit(row.names(mtx)[i], "_", fixed = T)
  source <- sapply(source, tail, 1)
  mtx.fac[i] <- source
}
mtx.fac <- as.factor(mtx.fac)
mtx.col <- colorRampPalette(brewer.pal(n = length(levels(mtx.fac)), name = "Set1"))(length(levels(mtx.fac)))[as.numeric(mtx.fac)]

# create 3D matrix
X = array(0, dim = c(10, 6))
for (i in 1 : length(levels(mtx.fac))) {
  cur_dim <- array(mtx[as.numeric(mtx.fac) == 2,], dim = c(10, 6))
  X <- abind(X, cur_dim, along = 3)
}

# plot dendrogram
col_palette <- colorRampPalette(brewer.pal(n = 9, name = "Reds"))
col_palette <- colorRampPalette(c("red", "white", "green"))

heatmap.2(mtx,
dendrogram = "none",
hclustfun = function(x) {hclust(x, method = "average")},
Colv = F,
Rowv = F,
margins = c(4, 16),
col = col_palette(15),
tracecol = NULL,
cellnote = mtx,
notecol = "black",
colRow = mtx.col,
colCol = "black"
)
