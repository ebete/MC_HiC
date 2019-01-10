#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

# constants
match <- 1
mismatch <- 1
gap <- 6
extend <- 1
error_rate <- 0.1

# functions
score_optimal <- function(x) {
  x * match
}
score_avg <- function(x) {
  x * ((1 - error_rate) * match -
    error_rate / 3 * mismatch -
    error_rate / 3 * gap * 2)
}
bt_scoring <- function(x) {
  20 + 8 * log(x)
}
bt2_scoring <- function(x) {
  20 + 8 * sqrt(x)
}

# apply functions
df <- data.frame(x = 1 : 500)
df$score_optimal <- score_optimal(df$x)
df$score_avg <- score_avg(df$x)
df$bt_scoring <- bt_scoring(df$x)
df$bt2_scoring <- bt2_scoring(df$x)

# plotting
df.melted <- melt(df, id.vars = "x")
score_plot <- ggplot(df.melted, aes(x = x, y = value, color = variable)) +
  geom_line() +
  xlab("Alignment length") +
  ylab("Alignment score") +
  ggtitle("Bowtie2 scoring functions") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(margin = margin(r = 0.5, l = 0.1, unit = "cm"))
  ) +
  guides(colour = guide_legend(nrow = 2))
score_plot
