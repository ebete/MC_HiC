#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
    library(ggplot2)
})

# constants
match <- 1
mismatch <- 1
gap <- 6
extend <- 1
error_rate <- 0.1

# functions
score_optimal <- function(x) { x * match}
score_avg <- function(x) { x * ((1 - error_rate) * match -
    error_rate / 3 * mismatch -
    error_rate / 3 * gap * 2)}
bt_scoring <- function(x) { 20 + 8 * log(x)}
bt2_scoring <- function(x) { 20 + 8 * sqrt(x)}

# plotting
ggplot(data.frame(x = 0), aes(x = x)) +
    xlab("Alignment length") +
    ylab("Alignment score") +
    ggtitle("Bowtie2 scoring functions") +
    theme_bw() +
    scale_color_brewer(palette = "Set1") +
    stat_function(fun = score_optimal, colour = "green") +
    stat_function(fun = score_avg, colour = "orange") +
    stat_function(fun = bt_scoring, colour = "purple") +
    stat_function(fun = bt2_scoring, colour = "blue") +
    xlim(1, 500)
