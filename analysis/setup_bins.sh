#!/bin/bash

module load bedtools2/2.27.0 samtools/1.9
cd /scratch/lmalins/paint/maize-by-michelangelo/data/maps
BINSIZE=1000000

# make genome file
grep "^@SQ" tinymap.sam | cut -f 2-3 | sed -e 's/SN://g' -e 's/LN://g' > genome_for_bedtools.txt

# make bins file
bedtools makewindows -g genome_for_bedtools.txt -w $BINSIZE > $BINSIZE'_win.bed'

# format of a BED-3 FILE
# chrome \t start \t end \n
# where start is 0-based coords

# convert sam to bam for using bedtools
samtools view -b -o tinymap.bam tinymap.sam

# count reads in non-overlapping bins
bedtools coverage -a $BINSIZE'_win.bed' -b tinymap.bam > $BINSIZE'_binned_read_counts.bed'