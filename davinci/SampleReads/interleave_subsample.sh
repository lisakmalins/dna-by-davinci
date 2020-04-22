#!/bin/bash
# Original: Kirk Amundson 15 March 2019
# Edit: Lisa Malins 28 June 2019

# Background:
# Reads downloaded from NCBI SRA SRR2960981 are ~100x coverage of the
# maize genome, 2x250 Illumina paired end reads.

# Aim:
# Count k-mers in a pseudo-random subset of SRR2960981.fastq representing ~34x coverage
# of the genome.

# LM 28 June 2019 edit:
# Redoing jellyfish count with higher coverage subsample.

# Steps:
# 1. Using the same pseudo-random number generator seed, subsample from forward and reverse mates
# 2. Interleave reads

# Reads downloaded with following commands:
# module load aspera-connect/3.5.1 sratoolkit/2.8.2-1
# prefetch --ascp-path '/software/aspera-connect/3.5.1/static/bin/ascp|/software/aspera-connect/3.5.1/static/etc/asperaweb_id_dsa.openssh' SRR2960981
# fastq-dump --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip /scratch/lmalins/sra/SRR2960981.sra
# Note: --split-3 option generates separate files for forward and reverse reads

# Step 1: Subsample from forward and reverse read mates to get ~34x coverage of the genome
# with paired-end reads. This is analogous to what Albert et al (2019) PNAS did.

# sample from forward mates using 85 as seed for pseudo-random number generator
module load seqtk/1.3
seqtk sample -s85 SRR2960981_pass_1.fastq 0.42 > 85seed_42xsub_SRR2960981-1.fastq

# sample from reverse mates using 85 as seed for pseudo-random number generator
seqtk sample -s85 SRR2960981_pass_2.fastq 0.42 > 85seed_42xsub_SRR2960981-2.fastq

# Step 2: Interleave reads. This is done using the shell script ./fastInterleaveFromUnzip.sh
./fastInterleaveFromUnzip.sh 85seed_42xsub_SRR2960981-1.fastq 85seed_42xsub_SRR2960981-2.fastq 85seed_42xsub_SRR2960981.fastq

# file 85seed_42xsub_SRR2960981.fastq can then be input for k-mer counting
# What is the required hash size for a set of 42x reads for 17mers and a 1% error rate?
# 2.5e9 + 2.5e9 * 42 * 0.01 * 17 = 2.035e10
# Roughly 20GB
# example command:
# jellyfish count 85seed_42xsub_SRR2960981.fastq -m 17 -C -s 20G -t 7 -o 85seed_42xsub_SRR2960981_17mer.jf

# If wanting to use a Bloom filter to keep memory usage down:
# Total k-mer expected more than once: 2.5e9
# Total k-mer in dataset: 2.035e10
# One pass method example command:
# jellyfish count 85seed_34xsub_SRR2960981.fastq.gz -m 17 -C -s 3G -t 7 -o 85seed_42xsub_SRR2960981_17mer.jf --bf-size 20G

# Two pass method example command:
# jellyfish bc 85seed_34xsub_SRR2960981.fastq -m 17 -C -s 20G -t 7 -o 85seed_42xsub_SRR2960981_17mer.bc
# jellyfish count 85seed_34xsub_SRR2960981.fastq.gz -m 17 -C -s 3G -t 7 -bc 85seed_42xsub_SRR2960981_17mer.bc -o 85seed_42xsub_SRR2960981_17mer.jf
