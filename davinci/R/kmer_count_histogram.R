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

# Find peak
# Ignore left tail
kmers_abridged <- kmers %>% slice(-1:-10)
# Store y-value of peak (greatest number of k-mers per abundance)
peak_y <- max(kmers_abridged$number)
# Store x-value of peak (abundance with greatest number of k-mers)
peak_x <- kmers_abridged[[which.max(kmers_abridged$number), "abundance"]]
print(paste("Peak detected at abundance", peak_x,
            "with", peak_y, "k-mers"))

# X-axis cutoff at 0.997 cumulative fraction
# Fixes x-axis stretched by tiny number of high-abundance k-mers
x_cutoff <- min(which(kmers$cum_fraction > 0.997))
# Y-axis cutoff will be 110% of y_peak
y_cutoff <- round(peak_y * 1.1)

# Plot
ggplot(filter(kmers, abundance >= 4), aes(x = abundance, y = number)) +
  geom_histogram(binwidth = 1, stat = "identity") +
  scale_x_continuous(limits = c(0, x_cutoff),
                     breaks = c(seq(from=0, to=x_cutoff, by=10),
                                peak_x)) +
  scale_y_continuous(limits = c(0, y_cutoff)) +
  geom_vline(xintercept = peak_x) +
  labs(title = "Number of K-mers by Abundance",
       x = "K-mer abundance in reads",
       y = "Number of distinct k-mers")


# Save
ggsave(output, plot = last_plot())
