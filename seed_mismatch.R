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
m <- 0 : 5
seed_len <- 30
mismatch_rate <- 0.1

# Create dataframe with probabilities
df <- data.frame(n = 1 : seed_len)
for (i in m) {
    tmp <- data.frame(dpm(i, s = seed_len, p = mismatch_rate))
    names(tmp) <- c(sprintf("mm_%d", i))
    df <- cbind(df, tmp)
}

# Plot probabilities
df <- melt(df, id.vars = "n")
ggplot(df, aes(x = n, y = value, colour = variable)) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1), labels = yscale) +
    ggtitle(sprintf("Seed discovery probability for ONT (error rate: %.0f%%)", mismatch_rate * 100)) +
    xlab("Seed length") +
    ylab("Discovery probability") +
    scale_colour_brewer(palette = "Set1")
