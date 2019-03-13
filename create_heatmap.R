#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
  library(gplots)
  library(RColorBrewer)
  library(abind)
})

# load data
df <- read.csv("/data0/thom/mc4c_aligntest2/filtered/fragments.csv", sep = "\t", header = T)
mtx <- data.matrix(df[, - 1])
row.names(mtx) <- df[, 1]
colnames(mtx) <- 1 : ncol(mtx)
# remove first few columns
merge_from <- 8
remove_to <- 2
merged_tail <- apply(mtx[, merge_from : ncol(mtx)], 1, sum)
mtx.subset <- cbind(mtx[, - (merge_from : ncol(mtx))], "8+" = merged_tail)
mtx.subset <- mtx.subset[, - (1 : remove_to)]

# sort based on total mapped fragments
mtx.scores <- matrix(mtx.subset[, 1]) * (remove_to + 1)
for (i in 2 : ncol(mtx.subset)) {
  mtx.scores <- mtx.scores + mtx.subset[, i] * (i + remove_to)
}
mtx.subset <- mtx.subset[order(mtx.scores[, 1], decreasing = T),]

# plot heatmap
col_palette <- colorRampPalette(brewer.pal(n = 3, name = "PuOr"))
heatmap.2(
mtx.subset,
  dendrogram = "none",
  Colv = F,
  Rowv = F,
key = T,
tracecol = NULL,
cellnote = mtx.subset,
notecol = "black",
  margins = c(4, 20),
  col = col_palette(15),
main = "Aligned fragments per read for different alignment strategies",
xlab = "Fragments per read"
)
