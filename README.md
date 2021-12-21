# DNA by da Vinci
Pipeline for designing whole-chromosome oligo-FISH paints. Inspired by Jiming Jiang and Jim Birchler's [chromosome-specific FISH library in maize](https://doi.org/10.1073/pnas.1813957116) (nicknamed "Maize by Monet").

Developed by Lisa Malins, supervised by Kirk Amundson in the [Luca Comai Lab](http://comailab.genomecenter.ucdavis.edu/index.php/Main_Page) at the [UC Davis Genome Center](https://genomecenter.ucdavis.edu/).

Check out my 3-minute intro to DNA by da Vinci and whole-chromosome paints on YouTube!  
[![Link to YouTube video, "DNA by da Vinci – Undergrad Slam"](https://img.youtube.com/vi/u5k5BFxsBy8/0.jpg)](https://www.youtube.com/watch?v=u5k5BFxsBy8)  
[DNA by da Vinci – Undergrad Slam - YouTube](https://www.youtube.com/watch?v=u5k5BFxsBy8)

## Summary
### __Part 1A__: Slice genome into long oligos
1. Slice a reference genome into overlapping long oligos: for example, all 45-mers separated by a step size of 3.
2. Map these oligos back to the reference genome to see where they would hybridize and how well.
3. Filter out oligos that map back to more than one location.
4. Filter out oligos that would make poor probes because they would not hybridize well or would form hairpins.

### __Part 1B__: Count short oligo frequencies
1. Given a library of sequencing reads, use Jellyfish to count the frequency of short oligos: for example, 17-mers.
2. Use Jellyfish to create a histogram of the k-mer frequencies. The major peak of this histogram represents the approximate coverage of the library.

### __Part 2__: Calculate scores to select probes
1. Calculate k-mer scores for each long oligo kept in the pipeline. For example, in each 45-mer, look up the counts for every 17-mer it contains. The k-mer score for the 45-mer will be the sum of the counts of its 17-mers.
2. Create a histogram of the k-mer scores for the long oligos.
3. Select the oligos with k-mer scores in the peak region of the histogram. These will be the best probes.

## Test Case
1. Git clone this repository.
```bash
git clone https://github.com/lisakmalins/dna-by-davinci.git
cd dna-by-davinci
```

2. Install Conda.
```bash
#TODO
```

3. Build and activate the Conda environment.
```bash
conda env create -f environment.yaml
source activate davinci
```

4. Go to data/seqs/ folder, download reference genome, unzip, and rename.
```bash
cd data/seqs
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
gunzip -c Zea_mays.AGPv4.dna.toplevel.fa.gz > Zea_mays.AGPv4.fa
```

5. Check the reference genome you downloaded for extraneous sequences you don't care about.

  For example, there may be plastid DNA or unincorporated contigs. Run the following command to see the headers for all your sequences.

  ```bash
grep "^>" Zea_mays.AGPv4.fa
```

  You should see 10 nuclear chromosomes, 2 plastids, and a bunch of contigs. Run the following command to save the headers for only the 10 chromosomes into `Zea_mays.AGPv4.keep`:

  ```bash
grep "^>" Zea_mays.AGPv4.fa | head -n 10 > Zea_mays.AGPv4.keep
```

  You could also use copy-paste with the command-line text editor Vim.

6. (Optional) If you have Aspera-connect on your machine, feel free to prefetch SRR2960981.sra and place in data/reads/ folder to make the fastq-dump step faster.

7. Run snakemake!
```bash
cd ../..  # Go back up to dna-by-davinci directory
snakemake --cores 16
```


## References
Albert PS, Zhang T, Semrau K, Rouillard JM, Kao YH, Wang CJ, Danilova TV, Jiang J, Birchler JA. Whole-chromosome paints in maize reveal rearrangements, nuclear domains, and chromosomal relationships. *Proceedings of the National Academy of Sciences*. 2019 Jan 29;116(5):1679-85.
