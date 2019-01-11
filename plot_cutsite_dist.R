#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
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
# calculate ECDF of DpnII cut site intervals
sites_ecdf <- data.frame(interval = seq(1, 500, by = 1))
for (chr in levels(sites_dist$chromosome)) {
  f <- ecdf(sites_dist[sites_dist$chromosome == chr,]$distance)
  df <- data.frame(f(sites_ecdf$interval))
  colnames(df) <- c(chr)
  sites_ecdf <- cbind(sites_ecdf, df)
}
rm(df, f, chr)

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
  theme_classic() +
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
  theme_classic() +
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

# combine plots
grid.arrange(int_plot, dens_plot, ncol = 1)

# plot ECDF of DpnII interval
sites_ecdf.melted <- melt(sites_ecdf, id.vars = c("interval"))
colnames(sites_ecdf.melted) <- c("interval", "chromosome", "ecdf")
ecdf_plot <- ggplot(sites_ecdf.melted, aes(x = interval, y = ecdf, colour = chromosome)) +
  geom_line() +
  geom_vline(xintercept = 256, linetype = "dotted", color = "red") + # theoretical average DpnII interval
  scale_x_continuous(limits = c(0, 500)) +
  scale_y_continuous(labels = function(x) { sprintf("%.0f%%", x * 100)}, limits = c(0, 1)) +
  ggtitle("ECDF of DpnII site intervals in mm9") +
  xlab("Interval (bases)") +
  ylab("Fraction of intervals shorter") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  annotate("text", x = 256, y = 0.4, label = "Theoretical average interval", colour = "red", vjust = - 1, hjust = 0, angle = 270)
ecdf_plot

# plot EDF of DpnII interval
interval_avg <- mean(sites_dist$distance[sites_dist$distance < 2000]) # (ignore 2k+ intervals)
edf_plot <- ggplot(sites_dist, aes(x = distance, fill = chromosome)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 256, linetype = "dashed", color = "darkred") + # theoretical DpnII interval
  geom_vline(xintercept = interval_avg, linetype = "dashed", color = "darkgreen") + # average DpnII interval
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous() +
  ggtitle("Distribution of DpnII site intervals in mm9") +
  xlab("Interval (bases)") +
  ylab("Density") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.ticks.y = element_blank(),
  axis.text.y = element_blank()
  ) +
  annotate("text", x = 256, y = 0.0025, label = "Theoretical interval: 256", colour = "darkred", vjust = - 1, hjust = 0.5, angle = 270) +
  annotate("text", x = interval_avg, y = 0.0025, label = sprintf("Average interval: %.0f", interval_avg), colour = "darkgreen", vjust = - 1, hjust = 0.5, angle = 270)
edf_plot
