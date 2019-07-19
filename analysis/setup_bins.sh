#!/bin/bash
# Purpose: Set up windows for binned probe counts
# Usage: bash setup_bins.sh {map input filename} {windows output filename} {binsize}
# Default binsize: 1000000 bases

module load bedtools2/2.27.0
MAP=$1
OUTPUT=$2
BINSIZE=${3:-1000000}

# Find end of headers
HEADEREND=`grep -m 1 -v "^@" -n $MAP | cut -f1 -d:`

# make genome file
head -n $HEADEREND $MAP | grep "^@SQ" | \
cut -f 2-3 | sed -e 's/SN://g' -e 's/LN://g' | \

# make bins file
bedtools makewindows -g /dev/fd/0 -w $BINSIZE > $2

# format of a BED-3 FILE
# chrom \t start \t end \n
# where start is 0-based coords
