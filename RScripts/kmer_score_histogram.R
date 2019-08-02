# Usage: rscript binned_coverage.R {input.txt} {output.jpg}
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]
score_histo <- read_csv(source)
ggplot(score_histo, aes(x = score, y = frequency)) +
  geom_histogram(binwidth = 1, stat="identity") +
  scale_x_continuous(limits = c(0,3000)) +
  labs(title = "K-mer Score Distribution", x = "K-mer score", y = "Frequency")
ggsave(output, plot = last_plot()) 