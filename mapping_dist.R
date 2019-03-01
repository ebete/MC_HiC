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
maplen_comparison <- read.delim("/data0/thom/mergemap_splitmap/mergemap.csv", header = T, stringsAsFactors = F)
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

mapdiff_plot <- ggplot(maplen_comparison, aes(x = effective_appended, fill = "Length difference")) +
  geom_histogram(binwidth = 5) +
#  scale_x_continuous(limits = c(NA, 1500)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  geom_vline(xintercept = median(maplen_comparison$effective_appended), color = "darkgreen") +
  annotate("text", x = median(maplen_comparison$effective_appended), y = 1500, label = sprintf("median: %.0f", mean(maplen_comparison$effective_appended)), colour = "darkgreen", vjust = - 1, hjust = 0, angle = 270) +
  ggtitle(sprintf("Difference of map length between %d alignments", nrow(maplen_comparison)))
mapdiff_plot

# plot mapping performance
perf <- data.frame(sample = NULL, size = NULL, improved = NULL, total = NULL)
for (f in Sys.glob("/data0/thom/splitmap/*.csv")) {
  df <- read.delim(f)
  fname <- strsplit(basename(f), ".", fixed = T)[[1]][1]
  sample <- strsplit(fname, "_", fixed = T)[[1]][2]
  extend <- strsplit(fname, "_", fixed = T)[[1]][4]
  perf <- rbind(perf, data.frame(
  sample = sample,
  size = extend,
  improved = sum(df$event == "improved"),
  total = nrow(df)
  ))
}
perf$size <- reorder(perf$size, as.numeric(perf$size))

map_perf <- ggplot(perf, aes(x = sample, y = improved / total, fill = size)) +
  geom_col(position = "dodge") +
  theme_classic2() +
  labs(
  title = "Second iteration mapping performance",
  subtitle = "Mapping performance of the merged fragments",
  x = "Dataset",
  y = "Fraction of cases"
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0, 0, 0.1)) +
  geom_text(
  aes(label = improved),
  vjust = 1,
  hjust = - 0,
  angle = 90,
  col = "black",
  fontface = "bold",
  position = position_dodge(width = 1)
  ) +
  theme(
  plot.title = element_text(face = "bold", hjust = 0.5),
  axis.text.x = element_text(angle = 0, vjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set1")
map_perf

# raw read stats
read_stats <- read_ods("/data0/thom/mc4c_fa/stats.ods")
read_stats$decrease <- 1 - read_stats$fragments / read_stats$`raw reads`
read_stats.m <- melt(read_stats, id.vars = c("identifier", "decrease"), measure.vars = c("raw reads", "fragments"))
readstat_plot <- ggplot(read_stats.m, aes(x = identifier, y = value, fill = variable)) +
  geom_col(position = "identity", alpha = 1) +
  theme_classic2() +
  scale_y_continuous(expand = c(0, 0, 0.1, 0), labels = scales::scientific) +
  labs(
  title = "Number of reads per MC-4C dataset",
  x = "Dataset",
  y = "Read count"
  ) +
  geom_text(
  aes(label = value),
  vjust = - 0.5,
  hjust = 0.5,
  col = "black"
  ) +
  geom_text(
  data = read_stats.m[read_stats.m$variable == "fragments",],
  aes(label = sprintf("-%.1f%%", decrease * 100)),
  vjust = 1.3,
  hjust = 0.5,
  col = "yellow",
  fontface = "bold"
  ) +
  theme(
  plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set1")
readstat_plot

# MergeMap MAPQ scores
mergemap <- read.delim("/data0/thom/mergemap_splitmap/mergemap.csv", stringsAsFactors=FALSE)
ggplot(mergemap, aes(x=mapq)) +
  geom_histogram(binwidth = 1, aes(y=..density..)) +
  scale_y_continuous(expand = c(0,0), labels=percent) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = sprintf("Mapping quality of %d MergeMap alignments", nrow(mergemap)),
    x = "MAPQ",
    y = "Fraction of alignments"
  )
