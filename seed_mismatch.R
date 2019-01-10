#!/usr/bin/env Rscript

# Package import
suppressPackageStartupMessages({
  library(scales)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
})

# Function declarations
yscale <- function(x) { sprintf("%.0f%%", x * 100)}
dpm <- function(x, s=20, p=0.1) { 1 - pbinom(x, size = 1 : s, prob = p, lower.tail = F)}

# Constants declarations
m <- seq.int(0, 10, 2)
seed_len <- 100
mismatch_rate <- 0.1
label_points <- c(25, 50, 100)

# Create dataframe with probabilities
df <- data.frame(n = 1 : seed_len)
for (i in m) {
    tmp <- data.frame(dpm(i, s = seed_len, p = mismatch_rate))
  names(tmp) <- c(sprintf("%d", i))
    df <- cbind(df, tmp)
}

# Plot probabilities
df <- melt(df, id.vars = "n")
colnames(df) <- c("length", "k", "probability")
mm_plot <- ggplot(df, aes(x = length, y = probability, colour = k)) +
  geom_line() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1), labels = yscale) +
  labs(
  title = "Probabilities for at most k mismatches",
  subtitle = sprintf("ONT error rate: %.0f%%", mismatch_rate * 100)
  ) +
  xlab("Sequence length") +
  ylab("Probability") +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "right",
#    legend.title = element_blank(),
    legend.text = element_text(margin = margin(r = 0.5, l = 0.1, unit = "cm"))
  ) +
  guides(colour = guide_legend(ncol = 1)) +
  geom_text(aes(label = ifelse(length %in% label_points, sprintf("%.0f%%", probability * 100), "")), hjust = 1, vjust = 1.5, show.legend = F) +
  geom_point(data = df[df$length %in% label_points,], aes(colour = k), show.legend = F)
mm_plot

#ggarrange(mm_plot, score_plot, ncol = 2, labels = c("(1)", "(2)"))
