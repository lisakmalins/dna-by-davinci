# Usage: Rscript kmer_score_histogram.R {input.txt} {output.ext}

library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer score histogram from", source))
print(paste("Saving k-mer score plot to", output))

# Read in data
score_histo <- read_csv(source)

# Plot
ggplot(score_histo, aes(x = score, y = frequency)) +
  geom_col() +
  scale_x_continuous(limits = c(0,3000)) +
  labs(title = "K-mer Score Distribution", x = "K-mer score", y = "Frequency")

# Save
ggsave(output, plot = last_plot())
