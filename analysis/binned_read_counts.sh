#!/bin/bash

module load bedtools2/2.27.0 samtools/1.9
cd /scratch/lmalins/paint/maize-by-michelangelo/data/maps
BINSIZE=1000000

# convert sam to bam for using bedtools
samtools view -b -o tinymap.bam tinymap.sam

# count reads in non-overlapping bins
bedtools coverage -a $BINSIZE'_win.bed' -b tinymap.bam > $BINSIZE'_binned_read_counts.bed'
