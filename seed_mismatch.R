#!/usr/bin/env Rscript

# Package import
suppressPackageStartupMessages({
  library(scales)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(latex2exp)
})

# Function declarations
yscale <- function(x) { sprintf("%.0f%%", x * 100)}
dpm <- function(n, k=20, p=0.1) { 1 - pbinom(n, size = 1 : k, prob = p, lower.tail = F)}
make_plot <- function(n, k, p, label_points) {
  # n: sequence length; k: max errors; p: error probability

  # Create dataframe with probabilities
  df <- data.frame(n = 1 : n)
  for (i in k) {
    tmp <- data.frame(dpm(i, n, p))
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
    title = "Probabilities for at most k errors",
    subtitle = sprintf("ONT error rate: %.0f%%", p * 100)
    ) +
    xlab("Sequence length") +
    ylab(TeX("$Pr(X \\leq k)$")) +
    scale_colour_brewer(palette = "Dark2") +
    theme_classic() +
    theme(
    plot.title = element_blank(), #element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_blank(), #element_text(hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(margin = margin(r = 0.5, l = 0.1, unit = "cm"))
    ) +
    guides(colour = guide_legend(nrow = 1)) +
    geom_text(aes(label = ifelse(length %in% label_points, sprintf("%.0f%%", probability * 100), "")), hjust = 1, vjust = 1.5, show.legend = F) +
    geom_point(data = df[df$length %in% label_points,], aes(colour = k), show.legend = F)
  return(mm_plot)
}

# create plots
read_plot <- make_plot(500, seq.int(0, 50, 10), 0.1, seq.int(0, 500, 100))
seed_plot <- make_plot(20, c(0), 0.1, seq.int(0, 20, 5))

combined <- ggarrange(read_plot, seed_plot, ncol = 2, labels = c("(1)", "(2)"), hjust = 0, align = "v", common.legend = T, legend = "bottom")
annotate_figure(combined, top = "Probability for at most k errors")
