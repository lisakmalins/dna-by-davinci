## GetOligos.py

This script takes a fasta genome assembly and divides it into oligos of a given length and separated by a given step size.

Accepts as arguments the source filename, k-mer size, step size, and output filename.

Output is in fasta format:
```
>[chromosome number]_[index]
[sequence]
```

### Usage:
```
python GetOligos.py {reference genome filename} {oligo size} {step size} {output filename}
```

Example command: get 45-mers separated by step size of 3 from _Agra cadabra_ reference genome and save as `agra_cadabra_45mers.fa`
```
python GetOligos.py agra_cadabra_genome.fa 45 3 agra_cadabra_45mers.fa
```

GetOligos.py is dependent on k-mer and header classes. The header class will need to be rewritten for different assemblies since header formats vary.


### Kmer.py

GetOligos.py uses the Kmer class, which remembers its sequence, origin, and index within the chromosome. The class also has flags for end of section and end of file to return to main.

Methods:
* StartNewKMer() resets attributes and starts new k-mer from next section in file.
* GetNextKMer() forgets first *s* letters of sequence and appends next *s* letters or returns to main if end of section or end of file is reached.

### ZmaysB73Header.py

GetOligos.py also requires a header class to determine which sequences in a genome assembly to get k-mers from and to parse information out of the headers. This class was written for the Zea mays ssp mays cv B73 Reference Genome from <https://www.maizegdb.org/assembly>.

Since header format varies between genome assemblies, a different version of this file will need to be written for another assembly with a different header format.

Necessary components for rewrite:
* id -- string that holds identification for a given sequence (i.e., chromosome number, whether it's from a plastid or mitochondrion, or contig identification)
* EvalHeader() -- returns True or False whether main should get k-mers out of this sequence or not.
