# Usage: Rscript kmer_count_histogram.R {input.txt} {output.ext}

library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer count histogram from", source))
print(paste("Saving k-mer count plot to", output))

# Read in data
kmers <- read_delim(source, delim=" ", col_names = F, col_types='ii')

# Find peak
kmers_abridged <- kmers %>% slice(-1:-10)
peak <- kmers_abridged[which.max(kmers_abridged$X2), 1] %>% as.numeric()

# Plot
ggplot(filter(kmers, X1 >= 4), aes(x = X1, y = X2)) +
  geom_histogram(binwidth = 1, stat = "identity") +
  scale_x_continuous(limits = c(0, 150), breaks = c(peak, 0, 50, 100, 150)) +
  scale_y_continuous(limits = c(0, 2e7)) +
  labs(title = "Number of K-mers per Frequency", x = "Frequency", y = "Number of K-mers") +
  geom_vline(xintercept = peak)

# Save
ggsave(output, plot = last_plot())
