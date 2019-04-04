#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
  library(gplots)
  library(RColorBrewer)
  library(abind)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
})

# load data
df <- read.csv("/data0/thom/mc4c_vs_mchic/fragments.csv", sep = "\t", header = T)
# df <- read.csv("/data0/thom/mc4c_aligntest2/filtered/fragments.csv", sep = "\t", header = T)
mtx <- data.matrix(df[, - 1])
row.names(mtx) <- df[, 1]
colnames(mtx) <- 1 : ncol(mtx)
# remove first few columns
merge_from <- 8
remove_to <- 2
merged_tail <- apply(mtx[, merge_from : ncol(mtx)], 1, sum)
mtx.subset <- cbind(mtx[, - (merge_from : ncol(mtx))], "8+" = merged_tail)
mtx.subset <- mtx.subset[, - (1 : remove_to)]

# get total mapped fragments
mtx.scores <- matrix(mtx[, remove_to + 1]) * (remove_to + 1)
for (i in (remove_to + 2) : ncol(mtx)) {
  mtx.scores <- mtx.scores + mtx[, i] * i
}
# mtx.subset <- mtx.subset[order(mtx.scores[, 1], decreasing = T),]

df <- cbind(file = rownames(mtx.subset), as.data.frame(mtx.subset), total = mtx.scores)
df.m <- melt(df, id.vars = c("file", "total"))
df.m$value <- as.numeric(df.m$value)

# plot heatmap
# col_palette <- colorRampPalette(brewer.pal(n = 3, name = "PuRd"))
# heatmap.2(
# mtx.subset,
#   dendrogram = "none",
#   Colv = F,
#   Rowv = F,
# key = T,
# tracecol = NULL,
# cellnote = mtx.subset,
# notecol = "black",
#   margins = c(4, 20),
#   col = col_palette(15),
# main = "Aligned fragments per read for different alignment strategies",
# xlab = "Fragments per read"
# )

heatmap_plot <- df.m %>%
  mutate(file = reorder(file, total)) %>%
  group_by(file) %>%
  mutate(relcount = value / sum(value)) %>%
  ggplot(aes(y = file, x = variable, fill = relcount)) +
  geom_tile(show.legend = F) +
# stat_identity(geom = "text", aes(label = value), color = "white", fontface = "bold") +
  stat_identity(geom = "text", aes(label = sprintf("%.1f k", value / 1e3)), color = "white", fontface = "bold") +
  scale_x_discrete(expand = c(0 , 0)) +
  scale_y_discrete(expand = c(0 , 0)) +
  labs(
  title = "",
  x = "Mapped fragments per read",
  y = "Dataset"
  ) +
  theme_pubr(border = T) +
  theme(
  plot.title = element_text(hjust = 0.5),
  axis.text.y = element_text(colour = c("darkblue", "darkgreen", "darkgreen", "darkblue"))
  ) +
  scale_fill_gradient2(low = "#e7e1ef", mid = "#c994c7", high = "#dd1c77")
# heatmap_plot

getPalette = colorRampPalette(brewer.pal(9, "PuRd"))
bars_plot <- df.m %>%
  select(file, total) %>%
  distinct() %>%
  mutate(file = reorder(file, total)) %>%
  ggplot(aes(x = file, y = total, fill = file)) +
  geom_col(position = "stack", color = "black", width = 1) +
  geom_label(aes(y = 0, label = sprintf("%.1f k", total / 1e3)), vjust = 0.5, hjust = - 0.1, fontface = "bold", color = "black", fill = "white", alpha = 0.7) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0), labels = function(x) { sprintf("%.0f k", x / 1e3)}) +
  labs(
  title = "",
  x = "",
  y = "Total mapped fragments"
  ) +
  coord_flip() +
  theme_pubr(legend = "none") +
  theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
  ) +
  scale_fill_manual(values = getPalette(nrow(df)))
# bars_plot

combined <- ggarrange(heatmap_plot, bars_plot, ncol = 2, legend = "none", align = "h", widths = c(3, 2))
annotate_figure(combined, top = text_grob("Aligned fragments per read in MC Hi-C", face = "bold", size = 18))
