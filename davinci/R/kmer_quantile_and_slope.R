# 13 December 2020
# Lisa Malins
# kmer_quantile_and_slope.R

# Reads k-mer histogram output from Jellyfish (number vs. abundance),
# and adds columns for slope, cumulative sum, and cumulative fraction.
# The cumulative fraction column is identical to quantiles
# of the original list of abundances.

# Usage: Rscript kmer_quantile_and_slope.R {input.txt} {output.tsv}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]

print(paste("Reading k-mer count histogram from", source))
print(paste("Saving modified k-mer count histogram to", output))

# Read in data
kmers <- read_delim(source,
                    delim = " ",
                    col_names = c("abundance", "number"),
                    col_types = 'ii')

# Add columns for cumulative sum and cumulative percent
kmers <- kmers %>%
  # Cumulative sum of number column
  mutate(cum_sum = cumsum(number)) %>%
  # Cumulative fraction of number column
  mutate(cum_fraction = cum_sum / sum(number))

# Approximate first derivative
kmers <- kmers %>%
  mutate(slope=c(NA, diff(number)))

# Set slope as NA on last row
# The last row is a catchall bucket for high abundances
# so slope is not relevant
kmers[nrow(kmers), "slope"] <- NA

# Write histogram + additional info to disk as TSV
write_tsv(kmers, output)
