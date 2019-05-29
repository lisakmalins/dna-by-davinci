## LoadKmerDict.py
Reads 17-mers and counts from Jellyfish dump file into a nested dictionary. It is designed to be run in interactive mode. After loading the dictionary, you can calculate the size of the dictionary in memory or go on to calculate k-mer scores with the next program.

From python shell:

```
# Load k-mer dictionary in to memory
dumpfile = "{filename.fa}"
exec(open("LoadKmerDict.py").read())

# Further commands:
# Print all entries in dictionary
counts
# Output number of entries in dictionary
num_entries
# Output size of dictionary in bytes
sys.getsizeof(counts)
# Query count for a particular sequence
counts['{first 6 letters}']['{first 12 letters}']['{all 17 letters}']
```

## CalcKmerScores.py
Calculates k-mer scores for 45-mers using nested dictionary of 17-mers and counts loaded into memory by `LoadKmerDict.py`.

This program is designed to run in interactive mode, it can't accept filenames as arguments. Instead, the source file is set with the variable oligofile (to read/write in fasta format) or samfile (to read/write in sam format).

From python shell:
```
# Load k-mer scores from Jellyfish dump
dumpfile = "{filename.fa}"
exec(open("LoadKmerDict.py").read())

# Read k-mers from fasta or sam file (choose one)
oligofile = "{filename.fa}"
samfile = "{filename.sam}"

# Specify output filename (if none provided, defaults to standard out)
outputfile = "{filename}"

# Calculate k-mer scores
exec(open("CalcKmerScores.py").read())
```

## Test files
* `dump100.fa` is a tiny Jellyfish dump file for testing. It is the first 100 lines of a real Jellyfish dump file of 17-mers from maize.
* `fakedump.fa` is an artificial Jellyfish dump file with 17-mers containing only contiguous A's and G's.
* `fake45mers.fa` is an artificial file of 45-mers containing only contiguous A's and G's. Every 17-mer in `fake45mers.fa` can be found in `fakedump.fa`.
