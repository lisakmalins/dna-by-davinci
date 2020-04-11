# Lisa Malins
# Snakefile for Maize by Michelangelo pipeline

"""
Preparatory steps:
Download genome, place in data/genome/ directory, and rename to match genome wildcard if necessary
(Optional) If you are downloading reads from NCBI and wish to prefetch the reads sra file before running fastq-dump, place the sra file in the data/reads/ directory
"""

configfile: "config.yaml"

rule targets:
    input:
        # Jellyfish arm
        "flags/jellyfish.done",

        # Oligos arm
        "flags/oligos.done",

        # Plots
        "flags/plots.done"


###--------------------- Download reads ---------------------###

wildcard_constraints:
    read=config["reads"],
    # Override snakemake default: p can be empty string
    p=".*"

# Prefer quick_fastq_dump if sra file is prefetched, but do regular fastq_dump otherwise
ruleorder: quick_fastq_dump > fastq_dump

# Slow version sans prefetch
rule fastq_dump:
    output:
        "data/reads/{read}.fastq.gz"
    shell:
        "fastq-dump --outdir data/reads --gzip \
        --skip-technical --readids --read-filter pass --dumpbase \
        --split-spot --clip {wildcards.read}"

# Fast version if sra file is prefetched
rule quick_fastq_dump:
    input:
        "data/reads/{read}.sra"
    output:
        "data/reads/{read}.fastq.gz"
    shell:
        "fastq-dump --outdir data/reads --gzip \
        --skip-technical --readids --read-filter pass --dumpbase \
        --split-spot --clip {input}"

###------------------- Estimate coverage --------------------###

ruleorder: estimate_bases > estimate_bases_gz

# Estimate number of bases in fastq reads.
# Calculate from read line length and number of lines.
rule estimate_bases:
    input:
        "data/reads/{read}.fastq"
    output:
        "data/reads/{read}_readlength.txt",
        "data/reads/{read}_numlines.txt"
    shell: """
        echo "Estimating read length"
        head -n 1000 {input} | paste - - - - | cut -f 2 | wc -L > {output[0]}
        wc -l {input} > {output[1]}
   """

# Estimate number of bases in gzipped fastq reads.
# Calculate from read line length and number of lines.
rule estimate_bases_gz:
    input:
        "data/reads/{read}.fastq.gz"
    output:
        "data/reads/{read}_readlength.txt",
        "data/reads/{read}_numlines.txt"
    threads:
        config["subsampling"]["threads"]
    shell: """
        echo "Estimating read length"
        gunzip -c {input} | head -n 1000 | paste - - - - | cut -f 2 | wc -L > {output[0]}
        unpigz -p {threads} -c {input} | wc -l {output[1]}
    """

checkpoint estimate_coverage:
    input:
        "data/reads/{read}_readlength.txt",
        "data/reads/{read}_numlines.txt"
    output:
        "data/reads/{read}_numbases.txt",
        "data/reads/{read}_approx_coverage.txt"
    run:
        # Load read length and number of lines in fastq from last step
        with open(input[0], 'r') as infile:
            readlength = int(infile.read().strip())
        with open(input[1], 'r') as infile:
            numlines = int(infile.read().strip())

        # Approximate number of bases by typical length * number of reads
        numbases = readlength * numlines / 4
        with open(output[0], 'w') as out:
            out.write(str(round(numbases)) + "\n")

        # Approximate coverage by number of bases / genome size
        approx_coverage = numbases / config["genome_size"]
        with open(output[1], 'w') as out:
            out.write(str(round(approx_coverage)) + "\n")


###--------------------- Subsample if necessary ---------------------###
# Use unzipped reads if available
ruleorder: uninterleave > uninterleave_gz

# Uninterleave reads into forward and reverse, then gzip
rule uninterleave:
    input:
        "data/reads/{read}.fastq"
    output:
        temp("data/reads/{read}-1.fastq.gz"),
        temp("data/reads/{read}-2.fastq.gz")
    threads:
        config["subsampling"]["threads"]
    shell:
        """
        cat {input} \
            | paste - - - - - - - - \
            | tee >(cut -f 1-4 | tr '\t' '\n' | pigz -p {threads} > {output[0]}) \
            | cut -f 5-8 | tr '\t' '\n' | pigz -p {threads} > {output[1]}
        """

# Uninterleave gzipped reads into forward and reverse, then gzip
rule uninterleave_gz:
    input:
        "data/reads/{read}.fastq.gz"
    output:
        temp("data/reads/{read}-1.fastq.gz"),
        temp("data/reads/{read}-2.fastq.gz")
    threads:
        config["subsampling"]["threads"]
    shell:
        """
        unpigz -p {threads} -c {input} \
            | paste - - - - - - - - \
            | tee >(cut -f 1-4 | tr '\t' '\n' | pigz -p {threads} > {output[0]}) \
            | cut -f 5-8 | tr '\t' '\n' | pigz -p {threads} > {output[1]}
        """

# Subsample gzipped reads to a lower coverage, then gzip
rule subsample:
    input:
        "data/reads/{read}-{n}.fastq.gz"
    wildcard_constraints:
        n="1|2"
    params:
        seed=85,
        coverage=float(config["max_coverage"]) / 100
    output:
        temp("data/reads/{p}{read}-{n}.fastq.gz")
    threads:
        config["subsampling"]["threads"]
    shell:
        """
        unpigz -p {threads} -c {input} | \
        seqtk sample -s{params.seed} /dev/fd/0 {params.coverage} | \
        pigz -p {threads} > {output}
        """

# Interleave subsampled gzipped reads back into one gzipped file.
rule interleave:
    input:
        "data/reads/{p}{read}-1.fastq.gz",
        "data/reads/{p}{read}-2.fastq.gz"
    wildcard_constraints:
        # In this rule only, p must be non-empty
        p=".+"
    output:
        "data/reads/{p}{read}.fastq.gz"
    threads:
        config["subsampling"]["threads"]
    shell:
        """
        paste <(unpigz -p {threads} -c {input[0]} | paste - - - -) \
        <(unpigz -p {threads} -c {input[1]} | paste - - - -) \
        | tr '\t' '\n' | pigz -p {threads} > {output}
        """

###--------------------- Count 17-mer frequencies with Jellyfish ---------------------###

# Calculate expected k-mers from formula
def expected_kmers(wildcards=False):
    sys.stderr.write("Calculating hash size for jellyfish\n")

    G = config["genome_size"]
    c = config["max_coverage"]
    e = config["error_rate"]
    k = config["mer_size"]

    # Replace c with actual coverage if less than max
    try:
        readsfile = "data/reads/{read}_approx_coverage.txt".format(read=config["reads"])
        with open(readsfile, 'r') as f:
            approx_coverage = int(f.read().strip())
            if approx_coverage < int(config["max_coverage"]):
                c = approx_coverage
    except FileNotFoundError:
        sys.stderr.write("expected_kmers function says: could not open {}\n".format(readsfile))
        sys.stderr.write("Using value from config file instead: {}\n".format(c))

    exp_kmers = round(G + G * c * e * k)

    # Convert to gigabase or megabase for readability
    if exp_kmers > 1000000000:
        exp_kmers = str(round(exp_kmers / 1000000000)) + "G"
    elif exp_kmers > 1000000:
        exp_kmers = str(round(exp_kmers / 1000000)) + "M"

    return str(exp_kmers)

# Print calculated hash size by running `snakemake print_hash`
rule print_hash:
    run:
        print(expected_kmers())

# Temporarily unzip reads for Jellyfish
rule temp_unzip:
    input:
        "data/reads/{p}{read}.fastq.gz"
    output:
        temp("data/reads/{p}{read}.fastq")
    threads:
        config["jellyfish"]["threads"]
    shell:
        "unpigz -p {threads} -c {input} > {output}"

rule count_pass1:
    input:
        "data/reads/{p}{read}.fastq"
    params:
        k=config["mer_size"],
        bcsize=expected_kmers
    output:
        temp("data/kmer-counts/{p}{read}.bc")
    threads:
        config["jellyfish"]["threads"]
    shell:
        "jellyfish bc -m {params.k} -C -s {params.bcsize} -t {threads} -o {output} {input}"

rule count_pass2:
    input:
        fastq="data/reads/{p}{read}.fastq",
        bc="data/kmer-counts/{p}{read}.bc"
    params:
        k=config["mer_size"],
        genomesize=str(config["genome_size"])
    output:
        "data/kmer-counts/{p}{read}_{k}mer_counts.jf"
    threads:
        config["jellyfish"]["threads"]
    shell:
        "jellyfish count -m {params.k} -C -s {params.genomesize} -t {threads} --bc {input.bc} -o {output} {input.fastq}"

rule jellyfish_dump:
    input:
        "data/kmer-counts/{p}{read}_{k}mer_counts.jf"
    output:
        "data/kmer-counts/{p}{read}_{k}mer_dumps.fa"
    shell:
        "jellyfish dump {input} > {output}"

rule jellyfish_histo:
    input:
        "data/kmer-counts/{p}{read}_{k}mer_counts.jf"
    output:
        "data/kmer-counts/{p}{read}_{k}mer_histo.txt"
    shell:
        "jellyfish histo {input} > {output}"

# If coverage is too high and subsampling is necessary,
# these functions will use a prefix to request jellyfish on subsampled reads.
# If coverage is manageable, they will request jellyfish on the original reads.
def get_jelly_histo(wildcards):
    with open(checkpoints.estimate_coverage.get(read=config["reads"]).output[1]) as f:
        if int(f.read().strip()) > int(config["max_coverage"]):
            prefix = "85seed_{}sub".format(config["max_coverage"]) # TODO put seed in config
        else:
            prefix = ""
    return "data/kmer-counts/{p}{read}_{k}mer_histo.txt".format(p=prefix, read=config["reads"], k=config["mer_size"])

def get_jelly_dump(wildcards):
    with open(checkpoints.estimate_coverage.get(read=config["reads"]).output[1]) as f:
        if int(f.read().strip()) > int(config["max_coverage"]):
            prefix = "85seed_{}sub".format(config["max_coverage"])
        else:
            prefix = ""
    return "data/kmer-counts/{p}{read}_{k}mer_dumps.fa".format(p=prefix, read=config["reads"], k=config["mer_size"])

def get_jelly_histo_plots(wildcards):
    with open(checkpoints.estimate_coverage.get(read=config["reads"]).output[1]) as f:
        if int(f.read().strip()) > int(config["max_coverage"]):
            prefix = "85seed_{}sub".format(config["max_coverage"])
        else:
            prefix = ""
    return expand("data/plots/{p}{read}_{k}mer_histo.{ext}", \
    p=prefix, read=config["reads"], k=config["mer_size"], ext=["png", "pdf"])

rule jellyfish_done:
    input:
        get_jelly_histo,
        get_jelly_dump
    output:
        touch("flags/jellyfish.done")

rule calculate_peak:
    input:
        get_jelly_histo
    output:
        "data/kmer-counts/limits.txt"
    shell:
        "Rscript RScripts/calculate_limits.R {input} {output}"

###----------------------------- Download genome -----------------------------###
# rule download_genome:
#     output:
#         "data/genome/{{genome}}.{}"
#     shell:
#         """
#         wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
#         gunzip {genome}.fa.gz
#         """

###--------- Slice genome into overlapping 45-mers, map, and filter ----------###

# Handle genome to be either .fa or .fasta
GENOME, FASTA_EXT = config["genome"].rsplit(".", 1)
assert FASTA_EXT.lower() in ["fasta", "fa"], \
"Sorry, I don't recognize genome file extension {}.\n \
Make sure you put your genome assembly in data/genome/ \n \
and list the filename in config.yaml.\n \
Here's the genome filename currently listed in config.yaml:\n {}\n \
".format(FASTA_EXT, config["genome"])

rule relevant_seqs:
    output:
        temp("data/genome/{genome}_seqs.txt")
    run:
        with open(output[0], 'w') as f:
            f.write("\n".join(list(map(lambda x: str(x), config["sequences"] + [""] ))))


rule get_oligos:
    input:
        "data/genome/{{genome}}.{}".format(FASTA_EXT),
        "data/genome/{genome}_seqs.txt"
    output:
        "data/oligos/{genome}_45mers.fasta" # TODO wildcard mer
    log:
        "data/oligos/{genome}_45mers.log" # TODO wildcard mer
    shell:
        "python GetOligos/GetOligos.py {input} 45 3 {output} {log}" # TODO put params in config.yaml

# one option: split oligo reads here, then spawn parallel map/filter rules on split

rule bwa_index:
    input:
        "data/genome/{{genome}}.{}".format(FASTA_EXT),
    output:
        "data/genome/{{genome}}.{}.amb".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.ann".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.bwt".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.pac".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.sa".format(FASTA_EXT)
    shell:
        "bwa index {input}"

rule map_oligos:
    input:
        "data/genome/{{genome}}.{}.amb".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.ann".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.bwt".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.pac".format(FASTA_EXT),
        "data/genome/{{genome}}.{}.sa".format(FASTA_EXT),
        genome="data/genome/{{genome}}.{}".format(FASTA_EXT),
        oligos="data/oligos/{genome}_45mers.fasta"
    output:
        "data/maps/{genome}_45mers_unfiltered.sam"
    threads:
        config["mapping"]["threads"]
    shell:
        "bwa mem -t {threads} {input.genome} {input.oligos} > {output}"

rule filter_scatter:
    input:
        "data/maps/{genome}_45mers_unfiltered.sam"
    output:
        temp(expand("data/maps/split/{genome}_45mers_unfiltered_{chr}.sam", chr=["header"] + config["sequences"] , genome="{genome}"))
    run:
        # Just in case: cast config sequences to strings
        config_seqs = list(map(lambda x: str(x), config["sequences"]))

        # Set up dictionary of output files
        outfiles = dict()
        outfiles["header"] = open(output[0], 'w')
        for f in output[1:]:
            for s in config_seqs:
                if f.rsplit(".sam", 1)[0].split("unfiltered_")[-1] == s:
                    outfiles[s] = open(f, 'w')


        with open(input[0], 'r') as source:
            # Write headers
            sys.stderr.write("Splitting sam file into chunks for filtering by BWA scores.\n")
            sys.stderr.write("Writing headers\n")
            line = source.readline()
            assert line[0] == "@", "Expected header line beginning with @, instead found\n{}".format(line)
            while line[0] == "@":
                outfiles["header"].write(line)
                line = source.readline()

            current = "placeholder"

            sys.stderr.write("Writing sequences\n")
            while line:
                if line.split("_", 1)[0] != current:
                    current = line.split("_", 1)[0]
                    assert current in config_seqs, "Unrecognized sequence {}, line was\n{}".format(current, line)
                    sys.stderr.write("Writing sequence {}\n".format(current))
                while line.split("_", 1)[0] == current:
                    outfiles[current].write(line)
                    line = source.readline()


rule filter_oligos:
    input:
        "data/maps/split/{genome}_45mers_unfiltered_{n}.sam"
    output:
        temp("data/maps/split/{genome}_45mers_filtered_{n}.sam")
    params:
        bwa_min_AS=45,
        bwa_max_XS=31,
        primer3_min_TM=37,
        primer3_max_HTM=35,
        primer3_min_diff_TM=10,
        write_rejected=False
    script:
        "FilterOligos/FilterOligos.py"

rule filter_gather:
    input:
        "data/maps/split/{genome}_45mers_unfiltered_header.sam",
        expand("data/maps/split/{{genome}}_45mers_filtered_{chr}.sam", chr=config["sequences"])
    params:
        logs=expand("data/maps/split/{genome}_45mers_filtered_{chr}.log", genome=GENOME, chr=config["sequences"])
    output:
        "data/maps/{genome}_45mers_filtered.sam",
        "data/maps/{genome}_45mers_filtered.log"
    shell:
        """
        # Concatenate output
        for f in {input}; do cat $f >> {output[0]}; done
        # Concatenate logs
        for f in {params.logs}; do cat $f >> {output[1]}; done
        """

rule oligos_done:
    input:
        expand("data/maps/{genome}_45mers_filtered.{ext}", genome=GENOME, ext=["sam", "log"])
    output:
        touch("flags/oligos.done")

# could put a checkpoint for split sam files here. Merge before doing calc scores and the rest

###------------------- Calculate k-mer scores for 45-mers --------------------###

rule calc_scores:
    input:
        dump=get_jelly_dump,
        map="data/maps/{genome}_45mers_filtered.sam"
    log:
        "data/scores/{genome}_45mers_scores.log"
    output:
        "data/scores/{genome}_45mers_scores.sam"
    shell:
        "python CalcScores/CalcKmerScores.py {input.dump} {input.map} {output}"

rule score_histogram:
    input:
        "data/scores/{genome}_45mers_scores.sam"
    output:
        "data/scores/{genome}_45mers_scores_histo.txt"
    shell:
        "python ScoresHisto/ScoresHistogram.py {input} {output}"

rule score_select:
    input:
        "data/scores/{genome}_45mers_scores.sam",
        "data/kmer-counts/limits.txt"
    output:
        "data/probes/{genome}_45mers_probes_selected.sam"
    log:
        "data/probes/{genome}_45mers_probes_selected.log"
    script:
        "SelectScores/SelectScores.py"


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
        probes="data/probes/{genome}_45mers_probes_selected.sam",
        bins="data/coverage/{{genome}}_45mers_{binsize}_bins.bed".format(binsize=config["binsize"])
    output:
        "data/coverage/{genome}_45mers_probes_coverage.bed"
    shell:
        "bash analysis/binned_read_counts.sh {input.probes} {input.bins} {output}"


###-------------------------------- R plots ---------------------------------###

rule binned_count_plot:
    input:
        "data/coverage/{genome}_45mers_probes_coverage.bed"
    output:
        "data/plots/{genome}_45mers_probes_coverage.{ext}"
    shell:
        "Rscript RScripts/binned_coverage.R {input} {output}"

rule kmer_count_plot:
    input:
        "data/kmer-counts/{p}{read}_{k}mer_histo.txt"
    output:
        "data/plots/{p}{read}_{k}mer_histo.{ext}"
    wildcard_constraints:
        read=config["reads"]
    shell:
        "Rscript RScripts/kmer_count_histogram.R {input} {output}"

rule kmer_score_plot:
    input:
        "data/scores/{genome}_45mers_scores_histo.txt"
    output:
        "data/plots/{genome}_45mers_scores_histo.{ext}"
    shell:
        "Rscript RScripts/kmer_score_histogram.R {input} {output}"

rule plots_done:
    input:
        # Binned coverage plot
        expand("data/plots/{genome}_45mers_probes_coverage.{ext}", \
        genome=GENOME, ext = ["png", "pdf"]),

        # K-mer count histogram plot
        get_jelly_histo_plots,

        # K-mer score histogram plot
        expand("data/plots/{genome}_45mers_scores_histo.{ext}", \
        genome=GENOME, ext = ["png", "pdf"])

    output:
        touch("flags/plots.done")
