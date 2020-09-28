## FilterFasta.py
Python program to filter oligos that contain homopolymers or fail Primer3 criteria from a FASTA file. Inspiration from primer3_filter.py in [Chorus2](https://github.com/zhangtaolab/Chorus2) by [zhangtaolab](https://github.com/zhangtaolab).

- __Homopolymer filtering__: exclude oligos that contain the same nucleotide repeated for a given number of bases. (Default 5, set with `--homopolymer-length`)
- __Primer3 filtering__:
  - __Melting temperature filtering__: exclude oligos with a melting temperature that is too low, which would prevent the oligo from hybridizing to its target stably. (Default 37, set with `--min-tm`)
  - __Hairpin melting temperature filtering__: exclude oligos with a hairpin melting temperature that is too high, which would mean the oligo would tend to form a hairpin. (Default 35, set with `--max-htm`)
  - __Difference melting temperature filtering__: exclude oligos with insufficient difference between melting temperature and hairpin melting temperature. (Default difference 10, set with `--min-dtm`)

### Usage
```
usage: FilterFasta.py [-h] -i OLIGOS [-o OUTPUT] [--min-tm MIN_TM]
                      [--max-htm MAX_HTM] [--min-dtm MIN_DTM]
                      [--homopolymer-length HOMOPOLYMER_LENGTH] [--verbose]

Filter oligos from fasta. Check for homopolymers and primer3 criteria.

optional arguments:
  -h, --help            show this help message and exit
  -i OLIGOS, --in OLIGOS
                        input fasta filename
  -o OUTPUT, --output OUTPUT
                        output filename (default: standard out)
  --min-tm MIN_TM       minimum melting temperature (default: 37)
  --max-htm MAX_HTM     maximum hairpin melting temperature (default: 35)
  --min-dtm MIN_DTM     minimum difference between melting temperature and
                        hairpin melting temperature (default: 10)
  --homopolymer-length HOMOPOLYMER_LENGTH
                        minimum length of homopolymer to filter out (default:
                        5)
  --verbose             print filtered records and reason for filtering to
                        standard error (default: do not print)
```

## FilterSam.py
Python program to filter oligos that fail mapping criteria from a SAM file. Inspiration from bwa.py in [Chorus2](https://github.com/zhangtaolab/Chorus2) by [zhangtaolab](https://github.com/zhangtaolab).

- __Alignment score filtering__: exclude oligos with less than a minimum alignment score. (Default 45, set with `--bwa-max-AS`)
- __Suboptimal alignment score filtering__: exclude oligos with a suboptimal alignment score (refers to the score of the second best alignment) greater than a maximum. This filter is designed to remove multi-mapping reads. (Default 31, set with `--bwa-max-XS`)
- __Primer3 filtering__: Disabled by default.
  - A previous version of this script also filtered by Primer3 criteria. That functionality has been moved to `FilterFasta.py` so it is disabled by default, but it can be enabled with `--enable-primer3-filter`.

### Usage
```
usage: FilterSam.py [-h] -i INPUT [-o OUTPUT] [--bwa-min-AS MIN_AS]
                    [--bwa-max-XS MAX_XS] [--enable-primer3-filter]
                    [--min-TM MIN_TM] [--max-HTM MAX_HTM]
                    [--min-diff-TM MIN_DIFF_TM] [--write-rejects]

Filter oligos from SAM file based on BWA mapping statistics.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --in INPUT  input SAM filename
  -o OUTPUT, --output OUTPUT
                        args.output filename (default: standard out)
  --bwa-min-AS MIN_AS   minimum BWA alignment score (AS:i:) required to keep
                        oligo; suggested same as probe length (default: 45)
  --bwa-max-XS MAX_XS   maximum BWA suboptimal alignment score (XS:i:)
                        required to keep oligo; suggested probe length * 70%
                        (default: 31)
  --enable-primer3-filter
                        enable filtering by primer3 criteria
  --min-TM MIN_TM       minimum melting temperature (default: 37)
  --max-HTM MAX_HTM     maximum hairpin melting temperature (default: 35)
  --min-diff-TM MIN_DIFF_TM
                        minimum difference between melting temperature and
                        hairpin melting temperature (default: 10)
  --write-rejects       write rejected oligos to separate output file
```
