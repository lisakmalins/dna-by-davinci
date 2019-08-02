# Usage: rscript binned_coverage.R {input.bed} {chr #} {output.jpg}
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
source = args[1]
num = strtoi(args[2])
output = args[3]
#source = snakemake@input[1]
#num= = snakemake@params[1]
#output = snakemake@output[1]
coverage <- read_delim(source, delim="\t", col_names=F, col_types=cols(X1 = col_factor())) %>% 
  rename(chr = X1, start = X2, end = X3, feat_overlap = X4) %>%
  rename(base_overlap = X5, length = X6, frac_overlap = X7) %>%
  mutate(start / 1000000) %>% rename(pos = "start/1e+06")
ggplot(filter(coverage, chr == num), aes(x=pos, y=feat_overlap)) +
  geom_histogram(binwidth = 1, stat = "identity") +
  labs(x = "Position on chromosome (MB)", y = "Probe Count", title = paste("Chromosome", num))
ggsave(output, plot = last_plot()) 