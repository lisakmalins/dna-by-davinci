# Usage: rscript binned_coverage.R {input.bed} {output.ext}

library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]
print(paste("Reading coverage data from", source))
print(paste("Saving coverage plot to", output))

# Load coverage data
coverage <- read_delim(source, delim="\t", col_names=F, col_types=cols(X1 = col_factor())) %>%
  rename(chr = X1, start = X2, end = X3, feat_overlap = X4) %>%
  rename(base_overlap = X5, length = X6, frac_overlap = X7) %>%
  mutate(start / 1000000) %>% rename(pos = "start/1e+06")

# Plot for all chromosmoes with facet wrap
ggplot(coverage, aes(x=pos, y=feat_overlap)) +
  geom_histogram(binwidth = 1, stat = "identity") +
  labs(x = "Position on chromosome (MB)", y = "Probe Count", title = "Probe Coverage by Chromosome") +
  facet_wrap(vars(chr), ncol = 2, strip.position = "r")

# save
ggsave(output, plot = last_plot())
