# Usage: Rscript kmer_score_histogram.R {input.txt} {output.ext}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer score histogram from", source))
print(paste("Saving k-mer score plot to", output))

# Read in data
score_histo <- read_csv(source,
                        col_types="ii")

x_cutoff = 3000

# Plot
ggplot(score_histo, aes(x = score, y = frequency)) +
  # Histogram has already been computed so use geom_col() to display
  geom_col() +
  # Zoom x-axis without clipping
  coord_cartesian(xlim=c(0, x_cutoff)) +
  # Titles
  labs(title = "K-mer Score Distribution", x = "K-mer score", y = "Frequency")

# Save
ggsave(output, plot = last_plot())
