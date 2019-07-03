READS=["85seed_42xsub_SRR2960981"]

rule targets:
    input:
        expand("data/kmer-counts/{read}_17mer_histo.txt", read=READS),
        expand("data/kmer-counts/{read}_17mer_dumps.fa", read=READS),
        "data/maps/zmays_AGPv4_map_filtered_42x_scores_histo.txt"

#TODO reformat with wildcards
rule score_histogram:
    input:
        "data/maps/zmays_AGPv4_map_filtered_42x_scores.sam"
    output:
        "data/maps/zmays_AGPv4_map_filtered_42x_scores_histo.txt"
    shell:
        "python3 JellyfishKmers/ScoresHistogram.py {input} {output}"

#TODO redo this step in C++ or python script mdoe
# rule calc_scores:

# rule binned_counts:
#     input:
#     output:
#     shell:
#         "grep -m 1 -v "^@" -n zmays_AGPv4_map.sam | cut -f1 -d:"
#
# rule make_bins:
#     input:
#     output:
#     shell:

rule dump:
    input:
        "data/kmer-counts/{read}_17mer_counts.jf"
    output:
        "data/kmer-counts/{read}_17mer_dumps.fa"
    shell:
        "jellyfish dump {input} > {output}"

rule histo:
    input:
        "data/kmer-counts/{read}_17mer_counts.jf"
    output:
        "data/kmer-counts/{read}_17mer_histo.txt"
    shell:
        "jellyfish histo {input} > {output}"

rule count_pass2:
    input:
        fastq="data/reads/{read}.fastq",
        bc="data/kmer-counts/{read}.bc"
    output:
        "data/kmer-counts/{read}_17mer_counts.jf"
    threads: 16
    shell:
        "jellyfish count -m 17 -C -s 3G -t 16 --bc {input.bc} -o {output} {input.fastq}"

rule count_pass1:
    input:
        "data/reads/{read}.fastq"
    output:
        "data/kmer-counts/{read}.bc"
    threads: 16
    shell:
        "jellyfish bc -m 17 -C -s 20G -t 16 -o {output} {input}"
