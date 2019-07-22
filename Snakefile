READS=["85seed_42xsub_SRR2960981"]

# Limits for k-mer score filtering
countsper = 45 - 17 + 1 # each 45-mer contains this many 17-mers (29)
coverage = 34
lower = round(countsper * coverage * 0.375)
upper = round(countsper * coverage * 1.8125)

rule targets:
    input:
        expand("data/kmer-counts/{read}_17mer_histo.txt", read=READS),
        expand("data/kmer-counts/{read}_17mer_dumps.fa", read=READS),
        "data/scores/zmays_AGPv4_map_filtered_42x_scores_histo.txt",
        expand("data/scores/zmays_AGPv4_map_filtered_42x_scores_KS_{l}_{u}.sam", l=lower, u=upper)


### Calculate k-mer scores for 45-mers ###

rule calc_scores:
    input:
        dump=expand("data/kmer-counts/{read}_17mer_dumps.fa", read=READS),
        map="data/maps/zmays_AGPv4_map_filtered.sam"
    output:
        "data/scores/zmays_AGPv4_map_filtered_42x_scores.sam"
    shell:
        "python3 CalcScores/CalcKmerScores.py {input.dump} {input.map} {output}"

rule score_histogram:
    input:
        "data/scores/zmays_AGPv4_map_filtered_42x_scores.sam"
    output:
        "data/scores/zmays_AGPv4_map_filtered_42x_scores_histo.txt"
    shell:
        "python3 JellyfishKmers/ScoresHistogram.py {input} {output}"

rule score_filter:
    input:
        "data/scores/zmays_AGPv4_map_filtered_42x_scores.sam"
    output:
        "data/scores/zmays_AGPv4_map_filtered_42x_scores_KS_{lower}_{upper}.sam"
    shell:
        "python3 FilterByScore/FilterByScore.py {input} {lower} {upper} {output}"


### Count 17-mer frequencies with Jellyfish ###

rule count_pass1:
    input:
        "data/reads/{read}.fastq"
    output:
        "data/kmer-counts/{read}.bc"
    threads: 16
    shell:
        "jellyfish bc -m 17 -C -s 20G -t 16 -o {output} {input}"

rule count_pass2:
    input:
        fastq="data/reads/{read}.fastq",
        bc="data/kmer-counts/{read}.bc"
    output:
        "data/kmer-counts/{read}_17mer_counts.jf"
    threads: 16
    shell:
        "jellyfish count -m 17 -C -s 3G -t 16 --bc {input.bc} -o {output} {input.fastq}"

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


### Generate coverage plots ###

rule make_bins:
    input:
        "data/maps/zmays_AGPv4_map.sam"
    output:
        "data/maps/zmays_AGPv4_map_1000000_win.bed"
    shell:
        "bash setup_bins.sh {input} {output} 1000000"

rule binned_counts:
