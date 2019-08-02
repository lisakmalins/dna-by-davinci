# Usage: rscript binned_coverage.R {input.txt} {output.jpg}
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]
kmers <- read_delim(source, delim=" ", col_names = F, col_types='ii')
ggplot(filter(kmers, X1 >= 4), aes(x = X1, y = X2)) +
  geom_histogram(binwidth = 1, stat = "identity") +
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 34, 50, 100, 150)) +
  scale_y_continuous(limits = c(0, 2e7)) +
  labs(title = "Number of K-mers per Frequency", x = "Frequency", y = "Number of K-mers") +
  geom_vline(xintercept = c(15,34))
ggsave(output, plot = last_plot()) 