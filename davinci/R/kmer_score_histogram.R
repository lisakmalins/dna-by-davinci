# Usage: Rscript kmer_score_histogram.R {input.txt} {output.ext}
# A comma-separated list of coordinates to draw vertical lines may be given as a 3rd argument.
# For example: Rscript kmer_score_histogram.R {input.txt} {output.ext} 200,1000

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
#source = snakemake@input[1]
#output = snakemake@output[1]

print(paste("Reading k-mer score histogram from", source))
print(paste("Saving k-mer score plot to", output))

# Optional: accept comma-separated-list of vertical lines to draw on histogram
if (length(args) >= 3) {
  # Verify argument contains no illegal characters
  if (!grepl("^[0-9,]*$", args[3])) {
    stop(paste("Optional third argument to this script accepts\n",
               "a list of vertical lines to draw on the histogram.\n",
               "Must be a comma-separated list of integers.\n",
               "For example: 100,200,300\n",
               "Unexpected input was:", vlines_argument))
  } else {
    print("Vertical line argument detected.")
    print(paste("Vertical lines will be drawn at the following x-coordinates:", args[3]))
  }

  # Split argument into vector of integers
  vlines <- args[3] %>% strsplit(",") %>% unlist() %>% as.integer()

} else {
  # If no argument provided, save empty vector instead
  vlines <- vector(mode = "integer")
}

# Read in data
score_histo <- read_csv(source,
                        col_types="ii")

x_cutoff = 3000

# Plot
ggplot(score_histo, aes(x = score, y = frequency)) +
  # Histogram has already been computed so use geom_col() to display
  geom_col() +
  # Zoom x-axis without clipping
  coord_cartesian(xlim=c(0, x_cutoff)) +
  # Titles
  labs(title = "K-mer Score Distribution", x = "K-mer score", y = "Frequency") +
  # Draw vertical lines on histogram if provided
  # Does nothing if vlines is empty vector
  geom_vline(xintercept=vlines, color="#181818") +
  annotate(geom = "label",
           label = vlines,
           x = vlines,
           y = max(score_histo$frequency)*1.05)

# Save
ggsave(output, plot = last_plot())
