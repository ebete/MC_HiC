#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(RColorBrewer)
})

undigested <- scan("/data0/thom/temp/LVR_HS5_NP_len.txt", numeric())
digested <- scan("/tmp/hs5.txt", numeric())
undigested_aligned <- scan("/data0/thom/temp/bwa_run_len.txt", numeric())

df <- rbind(
data.frame(source = "undigested", length = undigested),
data.frame(source = "digested", length = digested),
data.frame(source = "undigested.aligned", length = undigested_aligned)
)

ggplot(df, aes(x = length, fill = source)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(10, 1e4), labels = function(x) {sprintf("%.2f Kb", x / 1e3)}, trans = "log10", breaks = trans_breaks("log10", inv = function(x) { 10 ^ x}, n = 10)) +
  scale_y_continuous() +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold"),
  legend.position = "bottom"
  ) +
  xlab("Read length") +
  ylab("Density") +
  ggtitle("Distribution of read lengths")
