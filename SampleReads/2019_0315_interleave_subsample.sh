#!/bin/bash
# Kirk Amundson
# 15 March 2019

# Background:
# Reads downloaded from NCBI SRA SRR2960981.fastq.gz are ~100x coverage of the
# maize genome, 2x250 Illumina paired end reads. This set of reads is too large to read
# into memory at once on an AWS server with 32G RAM.

# Aim:
# Count k-mers in a pseudo-random subset of SRR2960981.fastq.gz representing ~30x coverage
# of the genome.

# Steps:
# 1. Uninterleave reads (currently interleaved and gzipped)
# 2. Using the same pseudo-random number generator seed, subsample from forward and reverse mates
# 3. Interleave reads

# Step 1: Uninterleave reads. This is done using the shell script ./fastUninterleaveUnzip.sh
# inspiration for this script courtesy of Daniel Standage: https://biowize.wordpress.com/2015/03/26/the-fastest-darn-fastq-decoupling-procedure-i-ever-done-seen/
cd ~/reads
./fastUninterleaveUnzip.sh SRR2960981.fastq.gz

# Step 2: Subsample from forward and reverse read mates to get ~34x coverage of the genome
# with paired-end reads. This is analogous to what Albert et al (2019) PNAS did.

# sample from forward mates using 85 as seed for pseudo-random number generatorn
seqtk sample -s85 SRR2960981-1.fastq.gz 0.34 > 85seed_34xsub_SRR2960981-2.fastq

# sample from reverse mates using 85 as seed for pseudo-random number generator
seqtk sample -s85 SRR2960981-2.fastq.gz 0.34 > 85seed_34xsub_SRR2960981-2.fastq

# Step 3: Interleave reads. This is done using the shell script ./fastInterleaveFromUnzip.sh
./fastInterleaveFromUnzip.sh 85seed_34xsub_SRR2960981-2.fastq 85seed_34xsub_SRR2960981-2.fastq 85seed_34xsub_SRR2960981.fastq.gz

# file 85seed_34xsub_SRR2960981.fastq.gz can then be input for k-mer counting
# What is the required hash size for a set of 34x reads for 17mers and a 1% error rate?
# 2.5e9 + 2.5e9 * 30 * 0.01 * 17 = 1.525e10
# Roughly 16G, manageable on UCD servers and on an AWS m4.2xlarge
# example command:
# gunzip -c 85seed_34xsub_SRR2960981.fastq.gz | jellyfish count /dev/fd/0 -m 17 -C -s 16G -t 7 -o 85seed_34xsub_SRR2960981_17mer.jf

# If wanting to use a Bloom filter to keep memory usage down:
# Total k-mer expected more than once: 2.5e9
# Total k-mer in dataset: 1.525e10
# example command:
# gunzip -c 85seed_34xsub_SRR2960981.fastq.gz | jellyfish count /dev/fd/0 -m 17 -C -s 16G -t 7 -o 85seed_34xsub_SRR2960981_17mer.jf --bf-size 3G
