#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
  library(ggpubr)
})

# draw mapped fragments on top of the original read
aln_coords <- read.csv("/data0/thom/conservative_aln/aln_coords.csv", sep = ";", header = T)
plot.new(); plot.window(ylim = c(1, nrow(aln_coords) + 1), xlim = c(min(aln_coords$start) - 1, max(aln_coords$end) + 1))
for (i in 1 : nrow(aln_coords)) {
  line_col <- hsv(h = if (aln_coords[i,]$ref == "chr7")2 / 3 else 1, s = min(aln_coords[i,]$mapq / 20, 1), v = 1, alpha = 1)
  segments(aln_coords[i,]$start, i, aln_coords[i,]$end, i, col = line_col, lty = 1, lwd = max(aln_coords[i,]$mapq / 10, 1))
}

# plot the distribution of the distance to DpnII sites
cutsite_dist <- rbind(
cbind(MAPQ = "1", read.csv("/data0/thom/conservative_aln/LVR_HS5_NP_lowQ_chimeric.csv", sep = "\t", header = T)),
cbind(MAPQ = "60", read.csv("/data0/thom/conservative_aln/LVR_HS5_NP-digested_chimeric.csv", sep = "\t", header = T))
)
cutsite_dist.melt <- melt(cutsite_dist, id.vars = "MAPQ", measure.vars = c("start_dist", "end_dist"))
cutsite_plot <- ggplot(cutsite_dist.melt, aes(x = abs(value), fill = MAPQ, color = MAPQ)) +
  geom_density(alpha = 0.5, position = "identity", adjust = 1) +
  scale_x_continuous(limits = c(NA, 200)) +
  scale_y_continuous(expand = c(0, 0)) +
#  coord_cartesian(ylim = c(0, 1e5)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  xlab("Distance to DpnII") +
  ggtitle("Distribution of distance between fragment ends and DpnII")
cutsite_plot

# plot the distribution of mapped fragment lengths
map_len <- rbind(
data.frame(MAPQ = "1", length = scan("/data0/thom/conservative_aln/LVR_HS5_NP_lowQ_chimeric_len.txt", numeric())),
data.frame(MAPQ = "60", length = scan("/data0/thom/conservative_aln/LVR_HS5_NP-digested_chimeric_len.txt", numeric()))
)
maplen_plot <- ggplot(map_len, aes(x = length, fill = MAPQ, color = MAPQ)) +
  geom_density(alpha = 0.5, adjust = 1) +
  scale_x_continuous(limits = c(0, 200)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  xlab("Length") +
  ggtitle("Mapped fragment length distribution")
maplen_plot

# plot the distribution of read fragment coverage
coverage_hq <- read.csv("/data0/thom/conservative_aln/LVR_HS5_NP-digested_chimeric_coverage.csv", sep = ";", header = T)
coverage_lq <- read.csv("/data0/thom/conservative_aln/LVR_HS5_NP_lowQ_chimeric_coverage.csv", sep = ";", header = T)
coverage <- rbind(
data.frame(MAPQ = "1", coverage = coverage_lq$coverage),
data.frame(MAPQ = "60", coverage = coverage_hq$coverage)
)
coverage_plot <- ggplot(coverage, aes(x = coverage, fill = MAPQ, color = MAPQ)) +
#  geom_density(alpha = 0.7, adjust = 1/2) +
  geom_histogram(binwidth = 0.01, alpha = 0.7, position = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(expand = c(0, 0), labels = scales::scientific) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  xlab("Fraction of read fragments mapped") +
  ggtitle("Distribution of reads mapped")
coverage_plot

# show cis/trans mapping MAPQ distribution
trans_mappings <- read.csv("/data0/thom/conservative_aln/cis_trans_qual.csv", sep = "\t", header = T, as.is = T)
trans_mappings$reference[trans_mappings$reference != "chr7"] <- "trans"
trans_mappings$reference[trans_mappings$reference == "chr7"] <- "cis"
trans_mappings$reference <- as.factor(trans_mappings$reference)
cistrans_plot <- ggplot(trans_mappings, aes(x = MAPQ, fill = reference, color = reference)) +
  geom_histogram(aes(y = ..count.. / sum(..count..)), alpha = .5, binwidth = 1, position = "identity") +
  stat_ecdf(pad = F, geom = "step", direction = "vh") +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  ggtitle("Distribution of MAPQ for cis- and trans-mapped fragments") +
  ylab("Fraction of mapped fragments")
cistrans_plot

# normalised alignment score
norm_aln <- read.delim("/data0/thom/conservative_aln/norm_aln.csv", header = T, sep = '\t')
norm_aln$norm_aln_score <- (norm_aln$norm_aln_score - min(norm_aln$norm_aln_score)) / (max(norm_aln$norm_aln_score) - min(norm_aln$norm_aln_score))
norm_aln$mapq <- factor(as.integer(norm_aln$mapq / 10) * 10, labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60"))
alnscore_plot <- ggplot(norm_aln, aes(x = norm_aln_score, fill = mapq)) +
  geom_density(alpha = .3) +
  scale_x_continuous() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "RdYlBu") +
  scale_color_brewer(palette = "RdYlBu") +
  theme_classic() +
  ggtitle("Distribution of normalised alignment scores") +
  xlab("Normalised alignment score") +
  ylab("Density")
alnscore_plot

# combine all plots
ggarrange(cutsite_plot, maplen_plot, coverage_plot, cistrans_plot, alnscore_plot, ncol = 2, nrow = 3, labels = c("(1)", "(2)", "(3)", "(4)", "(5)"))

# plot the length distribution of the two mapping approaches
maplen_comparison <- read.delim("/data0/thom/mc4c_fa_from_merging/mergemap_cut.csv", header = T, stringsAsFactors = F)
maplen_comparison$len_diff <- maplen_comparison$mergemap_length - maplen_comparison$original_length
maplen_comparison$mq_diff <- maplen_comparison$mergemap_mapq - maplen_comparison$original_mapq
maplen.melt <- melt(maplen_comparison[, c(2, 5)], measure.vars = c("original_length", "mergemap_length"))
maplen_plot <- ggplot(maplen.melt, aes(x = value, fill = variable)) +
  geom_histogram(alpha = .5, binwidth = 10, position = "identity", aes(y = ..count.. / sum(..count..))) +
  scale_x_continuous(limits = c(NA, 1500)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  geom_vline(xintercept = median(maplen_comparison$original_length), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = median(maplen_comparison$mergemap_length), color = "darkred", linetype = "dashed") +
  annotate("text", x = median(maplen_comparison$original_length), y = 0.03, label = sprintf("median: %.0f", median(maplen_comparison$original_length)), colour = "blue", vjust = - 1, hjust = 0, angle = 270) +
  annotate("text", x = median(maplen_comparison$mergemap_length), y = 0.03, label = sprintf("median: %.0f", median(maplen_comparison$mergemap_length)), colour = "darkred", vjust = - 1, hjust = 0, angle = 270)
maplen_plot

mapdiff_plot <- ggplot(maplen_comparison, aes(x = len_diff, fill = "Length difference")) +
  geom_histogram(binwidth = 30) +
  scale_x_continuous(limits = c(NA, 1500)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  geom_vline(xintercept = median(maplen_comparison$len_diff), color = "darkgreen") +
  annotate("text", x = median(maplen_comparison$len_diff), y = 1500, label = sprintf("median: %.0f", mean(maplen_comparison$len_diff)), colour = "darkgreen", vjust = - 1, hjust = 0, angle = 270) +
  ggtitle(sprintf("Difference of map length between %d alignments", nrow(maplen_comparison)))
mapdiff_plot
