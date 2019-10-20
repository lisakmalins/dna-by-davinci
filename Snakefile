# Lisa Malins
# Snakefile for Maize by Michelangelo pipeline

"""
Preparatory steps:
Download genome, place in data/seqs/ directory, and rename to match genome wildcard if necessary
(Optional) If you are downloading reads from NCBI and wish to prefetch the reads sra file before running fastq-dump, place the sra file in the data/reads/ directory
"""

configfile: "config.yaml"

# Constrain wildcards based on config file
def constrain(arg):
    if type(arg) == "list":
        return "|".join(arg)
    else:
        return(arg)
wildcard_constraints:
    read=constrain(config["reads"]),
    p=constrain(config["prefix"])

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
        # Jellyfish arm
        "flags/jellyfish.done",

        # Oligos arm
        "flags/oligos.done",

        # Coverage plots
        expand("data/plots/{genome}_45mers_{p}{read}_scores_{lower}_{upper}_coverage.{{ext}}".format( \
        genome=config["genome"], p=config["prefix"], read=config["reads"], lower=lower, upper=upper),
        ext = ["png", "pdf"]),

        # K-mer count histogram plot
        expand("data/plots/{p}{read}_17mer_histo.{{ext}}".format( \
        p=config["prefix"], read=config["reads"]), ext = ["png", "pdf"]),

        # K-mer score histogram plot
        expand("data/plots/{genome}_45mers_{p}{read}_scores_histo.{{ext}}".format( \
        genome=config["genome"], p=config["prefix"], read=config["reads"]), ext = ["png", "pdf"])

###--------------------- Download reads ---------------------###

# Prefer quick_fastq_dump if sra file is prefetched, but do regular fastq_dump otherwise
ruleorder: quick_fastq_dump > fastq_dump

# Slow version sans prefetch
rule fastq_dump:
    output:
        "data/reads/{read}.fastq"
    shell:
        "fastq-dump --outdir data/reads \
        --skip-technical --readids --read-filter pass --dumpbase \
        --split-spot --clip {wildcards.read}"

# Fast version if sra file is prefetched
rule quick_fastq_dump:
    input:
        "data/reads/{read}.sra"
    output:
        "data/reads/{read}.fastq"
    shell:
        "fastq-dump --outdir data/reads \
        --skip-technical --readids --read-filter pass --dumpbase \
        --split-spot --clip data/reads/{wildcards.read}.sra"



###--------------------- Subsample if necessary ---------------------###
# Uninterleave reads into forward and reverse
rule uninterleave:
    input:
        "data/reads/{read}-1.fastq"
    output:
        temp("data/reads/{read}-1.fastq"),
        temp("data/reads/{read}-2.fastq")
    shell:
        "bash SampleReads/fastUninterleave.sh {input}"

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
rule get_oligos:
    input:
        "data/seqs/{genome}.fa"
    output:
        "data/oligos/{genome}_45mers.fa"
    shell:
        "python GetOligos/GetOligos.py {input} 45 3 {output}"

rule bwa_index:
    input:
        "data/seqs/{genome}.fa",
    output:
        "data/seqs/{genome}.fa.amb",
        "data/seqs/{genome}.fa.ann",
        "data/seqs/{genome}.fa.bwt",
        "data/seqs/{genome}.fa.pac",
        "data/seqs/{genome}.fa.sa"
    shell:
        "bwa index {input}"

rule map_oligos:
    input:
        "data/seqs/{genome}.fa.amb",
        "data/seqs/{genome}.fa.ann",
        "data/seqs/{genome}.fa.bwt",
        "data/seqs/{genome}.fa.pac",
        "data/seqs/{genome}.fa.sa",
        genome="data/seqs/{genome}.fa",
        oligos="data/oligos/{genome}_45mers.fa"
    output:
        "data/maps/{genome}_45mers_unfiltered.sam"
    threads:
        12
    shell:
        "bwa mem -t {threads} {input.genome} {input.oligos} > {output}"

rule filter_oligos:
    input:
        "data/maps/{genome}_45mers_unfiltered.sam"
    output:
        "data/maps/{genome}_45mers_filtered.sam"
    params:
        bwa_min_AS=45,
        bwa_max_XS=31,
        primer3_min_TM=37,
        primer3_max_HTM=35,
        primer3_min_diff_TM=10,
        write_rejected=False
    script:
        "FilterOligos/FilterOligos.py"


###------------------- Calculate k-mer scores for 45-mers --------------------###

rule calc_scores:
    input:
        dump="data/kmer-counts/{p}{read}_17mer_dumps.fa",
        map="data/maps/{genome}_45mers_filtered.sam"
    log:
        "data/scores/{genome}_45mers_{p}{read}_scores.log"
    output:
        "data/scores/{genome}_45mers_{p}{read}_scores.sam"
    shell:
        "python CalcScores/CalcKmerScores.py {input.dump} {input.map} {output}"

rule score_histogram:
    input:
        "data/scores/{genome}_45mers_{p}{read}_scores.sam"
    output:
        "data/scores/{genome}_45mers_{p}{read}_scores_histo.txt"
    shell:
        "python ScoresHisto/ScoresHistogram.py {input} {output}"

rule score_select:
    input:
        "data/scores/{genome}_45mers_{p}{read}_scores.sam"
    output:
        "data/scores/{genome}_45mers_{p}{read}_scores_{lower}_{upper}.sam"
    shell:
        "python SelectScores/SelectScores.py {input} {lower} {upper} {output}"


###-------------------- Generate coverage histogram data ---------------------###

# Make TSV file of all sequences and lengths in genome
rule make_windows1:
    input:
        "data/maps/{genome}_45mers_unfiltered.sam"
    output:
        temp("data/coverage/{genome}_allseqs.tsv")
    shell:
        """
        module load bedtools2/2.27.0

        # Find end of headers
        HEADEREND=`grep -m 1 -v "^@" -n {input} | cut -f1 -d:`

        # make genome file
        head -n $HEADEREND {input} | grep "^@SQ" | \
        cut -f 2-3 | sed -e 's/SN://g' -e 's/LN://g' > {output}
        """

# Select sequences listed in config file and discard others
rule make_windows2:
    input:
        "data/coverage/{genome}_allseqs.tsv"
    output:
        "data/coverage/{genome}_seqs.tsv"
    run:
        out = open(output[0], 'w')
        for line in open(input[0], 'r').readlines():
            sequences = list(map(lambda x: str(x), config["sequences"]))
            if line.split()[0] in sequences:
                out.write(line)

# Make bins according to binsize in selected sequences
rule make_windows3:
    input:
        "data/coverage/{genome}_seqs.tsv"
    params:
        binsize=config["binsize"]
    output:
        "data/coverage/{{genome}}_45mers_{binsize}_bins.bed".format(binsize=config["binsize"])
    shell:
        "bedtools makewindows -g {input} -w {params.binsize} > {output}"


rule binned_counts:
    input:
        map="data/scores/{genome}_45mers_{p}{read}_scores_{lower}_{upper}.sam",
        bins="data/coverage/{{genome}}_45mers_{binsize}_bins.bed".format(binsize=config["binsize"])
    output:
        "data/coverage/{genome}_45mers_{p}{read}_scores_{lower}_{upper}_coverage.bed"
    shell:
        "bash analysis/binned_read_counts.sh {input.map} {input.bins} {output}"

rule binned_count_plot:
    input:
        "data/coverage/{genome}_45mers_{p}{read}_scores_{lower}_{upper}_coverage.bed"
    output:
        "data/plots/{genome}_45mers_{p}{read}_scores_{lower}_{upper}_coverage.{ext}"
    shell:
        "Rscript RScripts/binned_coverage.R {input} {output}"

rule kmer_count_plot:
    input:
        "data/kmer-counts/{p}{read}_17mer_histo.txt"
    output:
        "data/plots/{p}{read}_17mer_histo.{ext}"
    shell:
        "Rscript RScripts/kmer_count_histogram.R {input} {output}"

rule kmer_score_plot:
    input:
        "data/scores/{genome}_45mers_{p}{read}_scores_histo.txt"
    output:
        "data/plots/{genome}_45mers_{p}{read}_scores_histo.{ext}"
    shell:
        "Rscript RScripts/kmer_score_histogram.R {input} {output}"
