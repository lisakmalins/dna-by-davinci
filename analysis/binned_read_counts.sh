#!/bin/bash
# Purpose: Count probes in non-overlapping bins
# Usage: bash binned_read_counts.sh {}

module load bedtools2/2.27.0 samtools/1.9
SAM=$1
WINDOWS=$2
OUTPUT=$3

# convert sam to bam for using bedtools
samtools view -b $SAM | \

# count reads in non-overlapping bins
bedtools coverage -a $WINDOWS -b /dev/fd/0 > $OUTPUT
