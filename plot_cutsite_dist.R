#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

# import and reorder data
sites <- read.csv("~/mc_hic/cut_sites.csv", sep = ";")
sites$chromosome <- factor(sites$chromosome, levels(sites$chromosome)[c(1, 12 : 19, 2 : 11, 21, 22, 20)])
sites_dist <- read.csv("~/mc_hic/cut_sites_dist.csv", sep = ";")
sites_dist$chromosome <- factor(sites_dist$chromosome, levels(sites_dist$chromosome)[c(1, 12 : 19, 2 : 11, 21, 22, 20)])
# y-axis highlighting
axis_face_style <- rep("black", 22)
axis_face_style[16] <- "red"  # chr7

# plot DpnII intervals
int_plot <- ggplot(sites_dist, aes(x = distance, y = chromosome)) +
  geom_bin2d(bins = 100) +
  scale_x_continuous(limits = c(0, 1000), position = "top", breaks = pretty_breaks(n = 6)) +
#  scale_y_discrete(limits = "chr7") +
  scale_y_discrete(limits = rev(levels(sites_dist$chromosome))) +
  scale_fill_continuous(type = "viridis") + #trans="log"
  ggtitle("Distance between DpnII sites in mm9") +
  xlab("Site interval (bases)") +
  ylab("Chromosome") +
  theme_minimal() +
  theme(
  plot.title = element_text(hjust = 0.5),
  legend.position = "right",
  axis.text.y = element_text(color = axis_face_style)
  )

# plot DpnII density on chromosome
hbb_region <- data.frame(xmin = 105e6, xmax = 115e6, ymin = "chr6", ymax = "chr8")
dens_plot <- ggplot(sites, aes(x = position, y = chromosome)) +
  geom_bin2d(binwidth = 1e5) + # 100kb bins
  scale_x_continuous(position = "top", labels = function(x) { sprintf("%.0f Mb", x / 1e6)}) +
#  scale_y_discrete(limits = "chr7") +
  scale_y_discrete(limits = rev(levels(sites$chromosome))) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("DpnII sites in mm9") +
  xlab("Position (Mb)") +
  ylab("Chromosome") +
  theme_minimal() +
  theme(
  plot.title = element_text(hjust = 0.5),
  legend.position = "right",
  axis.text.y = element_text(color = axis_face_style)
  ) +
  geom_rect(data = hbb_region,
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  color = "red",
  alpha = 0,
  inherit.aes = FALSE
  )

grid.arrange(int_plot, dens_plot, ncol = 1)
