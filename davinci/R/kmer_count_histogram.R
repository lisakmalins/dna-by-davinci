# Usage: Rscript kmer_count_histogram.R {input.txt} {output.ext}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer count histogram from", source))
print(paste("Saving k-mer count plot to", output))

source("davinci/R/validate_kmer_histo_colnames.R")

# Read in data
kmers <- read_delim(source,
                    delim = "\t",
                    col_names = TRUE,
                    col_types = 'iiidi')

# Find slope global maximum.
# Should be located at inflection point on left side of main k-mer peak
slope_global_max <- kmers %>%
  # Exclude last row of histogram,
  #   which is a "catchall" for super-duper high abundance k-mers
  # If there are lots, the slope will spike on the very last row
  slice(kmers %>% head(n=-1) %>% pull(slope) %>% which.max()) %>%
  transmute(x=abundance, y=slope) %>%
  as.list()

print(paste("Slope global maximum detected at x-coordinate",
            slope_global_max$x))
print(paste("Searching for most frequent k-mer abundance,",
            "excluding abundances",
            slope_global_max$x,
            "and below"))

# Find peak
peak <- kmers %>%
  # Ignore abundances to left of slope global maximum
  # (should be inflection point on left slope of hill)
  filter(abundance > slope_global_max$x) %>%
  # Grab the row with the highest number
  filter(number == max(number)) %>%
  # Should be only one, but just in case
  slice(1) %>%
  # Store x and y coordinates
  transmute(x=abundance, y=number) %>%
  as.list()

print(paste("Peak detected at abundance", peak$x,
            "with", peak$y, "k-mers"))

# Set x-axis cutoff at 4 times the peak x-coordinate
x_cutoff <- peak$x * 4
# Y-axis cutoff will be 110% of peak y-coordinate
y_cutoff <- round(peak$y * 1.1)

# Plot
ggplot(kmers, aes(x = abundance, y = number)) +
  # Use geom_col() because histogram is pre-calculated
  geom_col() +
  # Zoom x- and y- axes without clipping
  coord_cartesian(xlim = c(0, x_cutoff),
                  ylim = c(0, y_cutoff)) +
  # Add x ticks/labels by 10
  scale_x_continuous(breaks = seq(from = 0,
                                  to = x_cutoff,
                                  by = 10)) +
  # Add vertical line and label for peak
  geom_vline(xintercept = peak$x) +
  annotate(geom = "label",
           label = peak$x,
           x = peak$x,
           y = peak$y * 1.075) +
  labs(title = "Number of K-mers by Abundance",
       x = "K-mer abundance in reads",
       y = "Number of distinct k-mers")

# Save
ggsave(output, plot = last_plot())
