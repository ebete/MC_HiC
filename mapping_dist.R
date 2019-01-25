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
  segments(aln_coords[i,]$start, i, aln_coords[i,]$end, i, col = aln_coords[i,]$ref, lty = aln_coords[i,]$strand + 1, lwd = 2)
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

ggarrange(cutsite_plot, maplen_plot, coverage_plot, cistrans_plot, ncol = 2, nrow = 2, labels = c("(1)", "(2)", "(3)", "(4)"))
