# Lisa Malins
# Snakefile for Maize by Michelangelo pipeline

"""
Preparatory steps:
Download genome, place in data/seqs/ directory, and rename to match genome wildcard if necessary
(Optional) If you are downloading reads from NCBI and wish to prefetch the reads sra file before running fastq-dump, place the sra file in the data/reads/ directory
"""

# from os import path

READS=["SRR2960981"]
PREFIX=["85seed_42sub_"]
GENOMES=["Zea_mays.AGPv4"]

# To avoid infinite recursion: wildcard "read" is alphanumeric,
# wildcard "prefix" is alphanumeric plus an underscore
wildcard_constraints:
    read="[A-Za-z0-9]+",
    p="\d+seed_\d+sub_"

# Limits for k-mer score filtering
coverage = 34 #TODO write rule to calculate coverage instead of hardcoding
# How many 17-mers in one 45-mer?    45 - 17 + 1 = 29
# Ratios taken from Albert PS, et al. Whole-chromosome paints in maize reveal rearrangements, nuclear domains, and chromosomal relationships. Proceedings of the National Academy of Sciences. 2019 Jan 29;116(5):1679-85.
lower = round((45 - 17 + 1) * coverage * 0.375)
upper = round((45 - 17 + 1) * coverage * 1.8125)

rule targets:
    params:
        lb=lower,
        ub=upper
    input:
        expand("data/kmer-counts/{p}{read}_17mer_histo.txt", p=PREFIX, read=READS),
        expand("data/kmer-counts/{p}{read}_17mer_dumps.fa", p=PREFIX, read=READS),
        expand("data/scores/{genome}_45mers_{p}{read}_scores_histo.txt", zip, genome=GENOMES, p=PREFIX, read=READS),
        expand("data/coverage/{genome}_45mers_{p}{read}_scores_{lower}_{upper}_coverage.bed", zip, genome=GENOMES, p=PREFIX, read=READS, lower=lower, upper=upper)

###--------------------- Download reads ---------------------###

# Prefer quick_fastq_dump if sra file is prefetched, but do regular fastq_dump otherwise
ruleorder: quick_fastq_dump > fastq_dump

# Slow version sans prefetch
rule fastq_dump:
    output:
        "data/reads/{read}_pass_1.fastq",
        "data/reads/{read}_pass_2.fastq"
    shell:
        "fastq-dump --outdir data/reads --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.read}"

# Fast version if sra file is prefetched
rule quick_fastq_dump:
    input:
        "data/reads/{read}.sra"
    output:
        "data/reads/{read}_pass_1.fastq",
        "data/reads/{read}_pass_2.fastq"
    shell:
        "fastq-dump --outdir data/reads --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip data/reads/{wildcards.read}.sra"


###--------------------- Subsample if necessary ---------------------###
rule uninterleave:


rule subsample:
    input:
        "data/reads/{read}_pass_1.fastq",
        "data/reads/{read}_pass_2.fastq"
    params:
        seed=85,
        coverage=0.42
    output:
        "data/reads/{p}{read}_pass_1.fastq",
        "data/reads/{p}{read}_pass_2.fastq"
    shell:
        """
        seqtk sample -s{params.seed} {input[0]} {params.coverage} > {output[0]}
        seqtk sample -s{params.seed} {input[1]} {params.coverage} > {output[1]}
        """

rule interleave:
    input:
        "data/reads/{p}{read}_pass_1.fastq",
        "data/reads/{p}{read}_pass_2.fastq"
    output:
        "data/reads/{p}{read}.fastq"
    shell:
        "bash SampleReads/fastInterleaveFromUnzip.sh {input} {output}"

###--------------------- Count 17-mer frequencies with Jellyfish ---------------------###

rule count_pass1:
    input:
        "data/reads/{p}{read}.fastq"
    output:
        "data/kmer-counts/{p}{read}.bc"
    threads: 16
    shell:
        "jellyfish bc -m 17 -C -s 20G -t 16 -o {output} {input}"

rule count_pass2:
    input:
        fastq="data/reads/{p}{read}.fastq",
        bc="data/kmer-counts/{p}{read}.bc"
    output:
        "data/kmer-counts/{p}{read}_17mer_counts.jf"
    threads: 16
    shell:
        "jellyfish count -m 17 -C -s 3G -t 16 --bc {input.bc} -o {output} {input.fastq}"

rule jellyfish_dump:
    input:
        "data/kmer-counts/{p}{read}_17mer_counts.jf"
    output:
        "data/kmer-counts/{p}{read}_17mer_dumps.fa"
    shell:
        "jellyfish dump {input} > {output}"

rule jellyfish_histo:
    input:
        "data/kmer-counts/{p}{read}_17mer_counts.jf"
    output:
        "data/kmer-counts/{p}{read}_17mer_histo.txt"
    shell:
        "jellyfish histo {input} > {output}"

###----------------------------- Download genome -----------------------------###
# rule download_genome:
#     output:
#         "data/seqs/{genome}.fa"
#     shell:
#         """
#         wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
#         gunzip {genome}.fa.gz
#         """

###--------- Slice genome into overlapping 45-mers, map, and filter ----------###
#TODO remove header class dependency
rule get_oligos:
    input:
        "data/seqs/{genome}.fa"
    output:
        "data/oligos/{genome}_45mers.fa"
    shell:
        "python {input} 45 3 {output}"

rule map_oligos:
    input:
        genome="data/seqs/{genome}.fa",
        oligos="data/oligos/{genome}_45mers.fa"
    output:
        "data/maps/{genome}_45mers_unfiltered.sam"
    threads:
        12
    shell:
        "bwa mem -t {threads} {input.genome} {input.oligos} > output"

#TODO make the bwa and primer3 thresholds into parameters
rule filter_oligos:
    input:
        "data/maps/{genome}_45mers_unfiltered.sam"
    output:
        "data/maps/{genome}_45mers_filtered.sam"
    shell:
        "python FilterOligos/FilterOligos.py {input} {output}"


###------------------- Calculate k-mer scores for 45-mers --------------------###

rule calc_scores:
    input:
        dump="data/kmer-counts/{p}{read}_17mer_dumps.fa",
        map="data/maps/{genome}_45mers_filtered.sam"
    output:
        "data/scores/{genome}_45mers_{p}{read}_scores.sam"
    shell:
        "python3 CalcScores/CalcKmerScores.py {input.dump} {input.map} {output}"

rule score_histogram:
    input:
        "data/scores/{genome}_45mers_{p}{read}_scores.sam"
    output:
        "data/scores/{genome}_45mers_{p}{read}_scores_histo.txt"
    shell:
        "python3 ScoresHisto/ScoresHistogram.py {input} {output}"

rule score_select:
    input:
        "data/scores/{genome}_45mers_{p}{read}_scores.sam"
    output:
        "data/scores/{genome}_45mers_{p}{read}_scores_{lower}_{upper}.sam"
    shell:
        "python3 SelectScores/SelectScores.py {input} {lower} {upper} {output}"


###-------------------- Generate coverage histogram data ---------------------###

rule make_bins:
    input:
        "data/maps/{genome}_45mers_unfiltered.sam"
    params:
        binsize=1000000
    output:
        "data/coverage/{genome}_45mers_{params.binsize}_bins.bed"
    shell:
        "bash analysis/setup_bins.sh {input} {output} {params.binsize}"

rule binned_counts:
    input:
        map="data/scores/{genome}_45mers_{p}{read}_scores_{lower}_{upper}.sam",
        bins="data/coverage/{genome}_45mers_{params.binsize}_bins.bed"
    params:
        binsize=1000000
    output:
        "data/coverage/{genome}_45mers_{p}{read}_scores_{lower}_{upper}_coverage.bed"
    shell:
        "bash analysis/binned_read_counts.sh {input.map} {input.bins} {output}"
