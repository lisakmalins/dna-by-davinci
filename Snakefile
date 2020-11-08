# Lisa Malins
# Snakefile for DNA by da Vinci pipeline

"""
Preparatory steps:
Download genome, place in data/genome/ directory, and rename to match genome wildcard if necessary
(Optional) If you are downloading reads from NCBI and wish to prefetch the reads sra file before running fastq-dump, place the sra file in the data/reads/ directory
"""

configfile: "config.yaml"

def selected_oligos_by_bed():
    if config.get("bed_regions") is not None:
        return ["data/probes/{genome}_{o}mers_probes_selected_{region}.bam".format(\
            genome=config["genome"].rsplit(".f", 1)[0], \
            o=config["oligo_size"], \
            region=region.rsplit(".bed")[0]) \
        for region in config.get("bed_regions")]
    else:
        return []

rule targets:
    input:
        # Jellyfish arm
        "flags/jellyfish.done",

        # Oligos arm
        "flags/oligos.done",

        # Plots
        "flags/plots.done",

        selected_oligos_by_bed()


###--------------------- Download reads ---------------------###

wildcard_constraints:
    read=config["reads"],
    # Override snakemake default: p can be empty string
    p=".*"

# Prefer prefetched_fastq_dump if sra file is prefetched, but do regular fastq_dump otherwise
ruleorder: parallel_prefetched_fastq_dump > parallel_fastq_dump

# Slow version sans prefetch
rule parallel_fastq_dump:
    output:
        "data/reads/{read}.fastq.gz"
    params:
        # Alternative temporary directory
        tmpdir="data/reads/tmp/",
        # Native fastq-dump options
        options="--outdir data/reads --gzip --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip"
    threads:
        config["fastq_dump"]["threads"]
    shell:"""
        mkdir -p {params.tmpdir}
        parallel-fastq-dump --sra-id {wildcards.read} -t {threads} --tmpdir {params.tmpdir} {params.options}
        # Remove `_pass` from the filename
        mv data/reads/{wildcards.read}_pass.fastq.gz data/reads/{wildcards.read}.fastq.gz
        """

# Fast version if sra file is prefetched
rule parallel_prefetched_fastq_dump:
    input:
        "data/reads/{read}.sra"
    output:
        "data/reads/{read}.fastq.gz"
    params:
        # Alternative temporary directory
        tmpdir="data/reads/tmp/",
        # Native fastq-dump options
        options="--outdir data/reads --gzip --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip"
    threads:
        config["fastq_dump"]["threads"]
    shell:
        "parallel-fastq-dump --sra-id {input} -t {threads} --tmpdir {params.tmpdir} {params.options}"

###------------------- Estimate coverage --------------------###

ruleorder: estimate_bases > estimate_bases_gz

# Estimate number of bases in fastq reads.
# Calculate from read line length and number of lines.
rule estimate_bases:
    input:
        "data/reads/{read}.fastq"
    output:
        "data/reads/{read}_numbases.txt"
    shell: """
        echo "Estimating number of bases in {input}"
        cat {input} | paste - - - - | cut -f 2 | wc -c > {output}
   """

# Estimate number of bases in gzipped fastq reads.
# Calculate from read line length and number of lines.
rule estimate_bases_gz:
    input:
        "data/reads/{read}.fastq.gz"
    output:
        "data/reads/{read}_numbases.txt"
    threads:
        config["subsampling"]["threads"]
    shell: """
        echo "Estimating number of bases in {input}"
        unpigz -p {threads} -c {input} | paste - - - - | cut -f 2 | wc -c > {output}
    """

checkpoint estimate_coverage:
    input:
        "data/reads/{read}_numbases.txt"
    output:
        "data/reads/{read}_approx_coverage.txt"
    run:
        # Load number of bases in fastq from last step
        with open(input[0], 'r') as infile:
            numbases = int(infile.read().strip())

        # Approximate coverage by number of bases / genome size
        approx_coverage = numbases / config["genome_size"]
        with open(output[0], 'w') as out:
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
        seed=config["subsampling"]["seed"],
        coverage=float(config["subsampling"]["max_coverage"]) / 100
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

###--------------------- Count k-mer frequencies with Jellyfish ---------------------###

# Calculate expected k-mers from formula
def expected_kmers(wildcards=False):
    sys.stderr.write("Calculating hash size for jellyfish\n")

    G = config["genome_size"]
    c = config["subsampling"]["max_coverage"]
    e = config["error_rate"]
    k = config["kmer_size"]

    # Replace c with actual coverage if less than max
    try:
        readsfile = "data/reads/{read}_approx_coverage.txt".format(read=config["reads"])
        with open(readsfile, 'r') as f:
            approx_coverage = int(f.read().strip())
            if approx_coverage < int(config["subsampling"]["max_coverage"]):
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

# Unzip reads for Jellyfish
rule unzip_reads:
    input:
        "data/reads/{p}{read}.fastq.gz"
    output:
        "data/reads/{p}{read}.fastq"
    threads:
        config["jellyfish"]["threads"]
    shell:
        "unpigz -p {threads} -c {input} > {output}"

rule count_pass1:
    input:
        "data/reads/{p}{read}.fastq"
    params:
        k=config["kmer_size"],
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
        k=config["kmer_size"],
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
def prefix():
    with open(checkpoints.estimate_coverage.get(read=config["reads"]).output[1]) as f:
        if int(f.read().strip()) > int(config["subsampling"]["max_coverage"]):
            return "{}seed_{}sub".format(
                    config["subsampling"]["seed"],
                    config["subsampling"]["max_coverage"])
        else:
            return ""

def get_jelly_histo(wildcards):
    return "data/kmer-counts/{p}{read}_{k}mer_histo.txt".format(
    p=prefix(), read=config["reads"], k=config["kmer_size"])

def get_jelly_dump(wildcards):
    return "data/kmer-counts/{p}{read}_{k}mer_dumps.fa".format(
    p=prefix(), read=config["reads"], k=config["kmer_size"])

def get_jelly_histo_plots(wildcards):
    return expand("data/plots/{p}{read}_{k}mer_histo.{ext}", \
    p=prefix(), read=config["reads"], k=config["kmer_size"], ext=["png", "pdf"])

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
        "Rscript davinci/R/calculate_limits.R {input} {output}"

###----------------------------- Download genome -----------------------------###
rule download_genome:
    output:
        "data/genome/{}".format(config["genome"])
    params:
        url=config["genome_download"]
    shell:
        """
        # Download genome from url and save to desired filename
        wget {params.url} --output-document {output}
        # If genome is gzipped; move to extension .gz and then unzip it
        if $(gzip -tq {output}); then mv {output} {output}.gz && gunzip {output}.gz; fi
        """
#TODO check if downloaded genome is gzipped
###--------- Slice genome into overlapping oligos, map, and filter ----------###

# Handle genome to be either .fa or .fasta
GENOME, FASTA_EXT = config["genome"].rsplit(".", 1)
assert FASTA_EXT.lower() in ["fasta", "fa", "fna", "fas"], \
"Sorry, I don't recognize genome file extension {}.\n \
Make sure you put your genome assembly in data/genome/ \n \
and list the filename in config.yaml.\n \
Here's the genome filename currently listed in config.yaml:\n {}\n \
".format(FASTA_EXT, config["genome"])


rule get_oligos:
    input:
        "data/genome/{{genome}}.{}".format(FASTA_EXT),
    output:
        "data/oligos/{genome}_{o}mers_{chr}.fasta"
    log:
        "data/oligos/{genome}_{o}mers_{chr}.log"
    # Fix ambiguous wildcard by prohibiting {chr} to end with 'filtered'
    wildcard_constraints:
        chr=".*(?<!filtered)"
    params:
        oligo_size=config["oligo_size"],
        step_size=config["step_size"]
    shell:
        "python davinci/GetOligos/GetOligos.py -g {input} \
        -m {params.oligo_size} -s {params.step_size} -o {output} -l {log}\
        --sequences {wildcards.chr}"

rule primer3_homopolymer_filter:
    input:
        "data/oligos/{genome}_{o}mers_{chr}.fasta"
    output:
        "data/oligos/{genome}_{o}mers_{chr}_filtered.fasta"
    params:
        min_tm=37,
        max_htm=35,
        min_dtm=10,
        max_homopolymer=5
    shell:
        "python davinci/FilterOligos/FilterFasta.py -i {input} -o {output} \
        --min-tm {params.min_tm} --max-htm {params.max_htm} --min-dtm {params.min_dtm} \
        --homopolymer-length {params.max_homopolymer}"

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
        oligos="data/oligos/{genome}_{o}mers_{chr}_filtered.fasta"
    output:
        "data/maps/split/{genome}_{o}mers_unfiltered_{chr}.sam"
    threads:
        config["mapping"]["threads"]
    shell:
        "bwa mem -t {threads} {input.genome} {input.oligos} > {output}"

rule filter_bwa:
    input:
        "data/maps/split/{genome}_{o}mers_unfiltered_{chr}.sam"
    output:
        "data/maps/split/{genome}_{o}mers_filtered_{chr}.sam"
    params:
        bwa_min_AS=45,
        bwa_max_XS=31
    shell:
        "python davinci/FilterOligos/FilterSam.py -i {input} -o {output} \
        --bwa-min-AS {params.bwa_min_AS} --bwa-min-XS {params.bwa_max_XS}"

rule merge_filtered:
    input:
        expand("data/maps/split/{{genome}}_{{o}}mers_filtered_{chr}.sam",
        chr=config["sequences"])
    output:
        "data/maps/{genome}_{o}mers_filtered.sam"
    shell:"""
        # Get @SQ headers only from first file
        grep "^@SQ" {input[0]} > {output}
        # Get non-header lines from all files
        for f in {input}; do grep -v "^@" $f >> {output}; done
        """
# Merging filtered oligo files is necessary while CalcScores.py reads scores into memory.
# If CalcScores uses a database instead, it could run in parallel.

rule oligos_done:
    input:
        expand("data/maps/{genome}_{o}mers_filtered.sam",
        genome=GENOME,
        o=config["oligo_size"])
    output:
        touch("flags/oligos.done")

# could put a checkpoint for split sam files here. Merge before doing calc scores and the rest

###------------------- Calculate k-mer scores for oligos --------------------###

rule calc_scores:
    input:
        dump=get_jelly_dump,
        map="data/maps/{genome}_{o}mers_filtered.sam"
    log:
        "data/scores/{genome}_{o}mers_scores.log"
    output:
        "data/scores/{genome}_{o}mers_scores.sam"
    shell:
        "python davinci/CalcScores/CalcKmerScores.py {input.dump} {input.map} {output}"

rule score_histogram:
    input:
        "data/scores/{genome}_{o}mers_scores.sam"
    output:
        "data/scores/{genome}_{o}mers_scores_histo.txt"
    shell:
        "python davinci/ScoresHisto/ScoresHistogram.py {input} {output}"

rule score_select:
    input:
        "data/scores/{genome}_{o}mers_scores.sam",
        "data/kmer-counts/limits.txt"
    output:
        "data/probes/{genome}_{o}mers_probes_selected.sam"
    log:
        "data/probes/{genome}_{o}mers_probes_selected.log"
    script:
        "davinci/SelectScores/SelectScores.py"


###-------------------- Generate coverage histogram data ---------------------###

# Make TSV file of all sequences and lengths in genome
rule make_windows1:
    input:
        "data/maps/{{genome}}_{o}mers_filtered.sam".format(o=config["oligo_size"])
    output:
        "data/coverage/{genome}_allseqs.tsv"
    shell:
        """
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
        "data/coverage/{{genome}}_{{o}}mers_{binsize}_bins.bed".format(
        binsize=config["binsize"])
    shell:
        "bedtools makewindows -g {input} -w {params.binsize} > {output}"


rule binned_counts:
    input:
        probes="data/probes/{genome}_{o}mers_probes_selected.sam",
        bins="data/coverage/{{genome}}_{{o}}mers_{binsize}_bins.bed".format(binsize=config["binsize"])
    output:
        "data/coverage/{genome}_{o}mers_probes_coverage.bed"
    shell:
        "bash davinci/BinnedCounts/binned_read_counts.sh {input.probes} {input.bins} {output}"


###-------------------------------- R plots ---------------------------------###

rule binned_count_plot:
    input:
        "data/coverage/{genome}_{o}mers_probes_coverage.bed"
    output:
        "data/plots/{genome}_{o}mers_probes_coverage.{ext}"
    shell:
        "Rscript davinci/R/binned_coverage.R {input} {output}"

rule kmer_count_plot:
    input:
        "data/kmer-counts/{p}{read}_{k}mer_histo.txt"
    output:
        "data/plots/{p}{read}_{k}mer_histo.{ext}"
    wildcard_constraints:
        read=config["reads"]
    shell:
        "Rscript davinci/R/kmer_count_histogram.R {input} {output}"

rule kmer_score_plot:
    input:
        "data/scores/{genome}_{o}mers_scores_histo.txt"
    output:
        "data/plots/{genome}_{o}mers_scores_histo.{ext}"
    shell:
        "Rscript davinci/R/kmer_score_histogram.R {input} {output}"

rule plots_done:
    input:
        # Binned coverage plot
        expand("data/plots/{genome}_{o}mers_probes_coverage.{ext}", \
        genome=GENOME, o=config["oligo_size"], ext = ["png", "pdf"]),

        # K-mer count histogram plot
        get_jelly_histo_plots,

        # K-mer score histogram plot
        expand("data/plots/{genome}_{o}mers_scores_histo.{ext}", \
        genome=GENOME,
        o=config["oligo_size"],
        ext = ["png", "pdf"])

    output:
        touch("flags/plots.done")


rule bedtools_intersect:
    input:
        sam="data/probes/{genome}_{o}mers_probes_selected.sam",
        bed="{region}.bed"
    output:
        "data/probes/{genome}_{o}mers_probes_selected_{region}.bam"
    shell:
        "samtools view -b {input.sam} | bedtools intersect -a - -b {input.bed} > {output}"
