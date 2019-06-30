#!/bin/bash
# fastInterleaveFromZip.sh
# 15 March 2019
# Interleave a pair of uninterleaved and gzipped fastq files, then gzip interleaved
# output is interleaved .fastq.gz

# USAGE:
# fastInterleave.sh reads-1.fastq reads-2.fastq reads.fastq

echo $1
echo $2
echo $3

paste <(zcat $1 | paste - - - -) <(zcat $2 | paste - - - -) | tr '\t' '\n' | gzip > $3
