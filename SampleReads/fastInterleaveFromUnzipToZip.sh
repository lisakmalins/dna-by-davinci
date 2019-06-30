#!/bin/bash
# fastInterleaveFromUnzipToZip.sh
# 15 March 2019
# Interleave a pair of uninterleaved and unzipped fastq files, then gzip interleaved
# output is interleaved .fastq.gz

# USAGE:
# fastInterleave.sh reads-1.fastq reads-2.fastq reads.fastq.gz

echo $1
echo $2
echo $3

paste <(cat $1 | paste - - - -) <(cat $2 | paste - - - -) | tr '\t' '\n' | gzip > $3
