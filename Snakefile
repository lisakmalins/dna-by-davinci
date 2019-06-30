READS=["85seed_42xsub_SRR2960981"]

rule all:
    input:
        expand("data/kmer-counts/{read}_17mer_histo.txt", read=READS),
        expand("data/kmer-counts/{read}_dumps.fa", read=READS)

rule dump:
    input:
        "data/kmer-counts/{read}_17mer_counts.jf"
    output:
        "data/kmer-counts/{read}_dumps.fa"
    shell:
        "jellyfish dump {input} > {output}"

rule histo:
    input:
        "data/kmer-counts/{read}_17mer_counts.jf"
    output:
        "data/kmer-counts/{read}_17mer_histo.txt"
    shell:
        "jellyfish histo {input}"

rule count_pass2:
    input:
        fastq="data/reads/{read}.fastq",
        bc="data/kmer-counts/{read}.bc"
    output:
        "data/kmer-counts/{read}_17mer_counts.jf"
    threads: 16
    shell:
        "jellyfish count -m 17 -C -s 3G -t 16 -bc {input.bc} -o {output} {input.fastq}"

rule count_pass1:
    input:
        "data/reads/{read}.fastq"
    output:
        "data/kmer-counts/{read}.bc"
    threads: 16
    shell:
        "jellyfish bc -m 17 -C -s 20G -t 16 -o {output} {input}"
