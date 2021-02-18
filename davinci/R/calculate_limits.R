# Usage: Rscript calculate_limits.R {input.txt} {output.tsv}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
if (length(args) >= 3) {
  peakoutput = args[3]
}
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer count histogram from", source))
print(paste("Saving k-mer count peak and k-mer score limits to", output))

source("davinci/R/validate_kmer_histo_colnames.R")

# Read in data
kmers <- read_delim(source,
                    delim="\t",
                    col_names = TRUE,
                    col_types='iiidi')

# Find k-mer count peak
kmers_abridged <- kmers %>% slice(-1:-10)
count_peak <- kmers_abridged[which.max(kmers_abridged$number), 1] %>% as.numeric()

# Decide score upper and lower bound
predicted_score_peak <- count_peak * (45 - 17 + 1)
score_lower_limit <- round(predicted_score_peak * 0.375)
score_upper_limit <- round(predicted_score_peak * 1.8125)

# Build matrix from key-value pairs and convert to dataframe
info <- rbind(
    c("count_peak", count_peak),
    c("predicted_score_peak", predicted_score_peak),
    c("score_lower_limit", score_lower_limit),
    c("score_upper_limit", score_upper_limit)
  ) %>% as.data.frame()

print(info)

# Write info to disk as 2-column TSV of key-value pairs
write_tsv(info, output, col_names=FALSE)
