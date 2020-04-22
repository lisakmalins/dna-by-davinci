# Usage: Rscript calculate_limits.R {input.txt} {output.txt}

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
if (length(args) >= 3) {
  peakoutput = args[3]
}
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer count histogram from", source))
print(paste("Saving k-mer count peak to", output))

# Read in data
kmers <- read_delim(source, delim=" ", col_names = F, col_types='ii')

# Find peak
kmers_abridged <- kmers %>% slice(-1:-10)
peak <- kmers_abridged[which.max(kmers_abridged$X2), 1] %>% as.numeric()
lower <- round((45 - 17 + 1) * peak * 0.375)
upper <- round((45 - 17 + 1) * peak * 1.8125)

# Write peak to disk
write(paste(peak, lower, upper), output)
