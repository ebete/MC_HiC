#!/usr/bin/env Rscript

# Package import
suppressPackageStartupMessages({
    library(scales)
    library(ggplot2)
    library(reshape2)
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
colnames(df) <- c("length", "mismatches", "probability")
ggplot(df, aes(x = length, y = probability, colour = mismatches)) +
  geom_line() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1), labels = yscale) +
  labs(
  title = "Probability of finding at most n mismatches",
  subtitle = sprintf("ONT error rate: %.0f%%", mismatch_rate * 100)
  ) +
  xlab("Sequence length") +
  ylab("Probability") +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  plot.subtitle = element_text(hjust = 0),
  legend.position = "right"
  ) +
  geom_text(aes(label = ifelse(length %in% label_points, sprintf("%.0f%%", probability * 100), "")), hjust = - 0.2, vjust = 0, show.legend = F) +
  geom_point(data = df[df$length %in% label_points,], aes(colour = mismatches), show.legend = F)
