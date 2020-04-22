#!/bin/bash
# Purpose: Count probes in non-overlapping bins
# Usage: bash binned_read_counts.sh {sam input file} {windows input file} {counts output file}

SAM=$1
WINDOWS=$2
OUTPUT=$3

# convert sam to bam for using bedtools
samtools view -b $SAM | \

# count reads in non-overlapping bins
bedtools coverage -a $WINDOWS -b /dev/fd/0 > $OUTPUT
