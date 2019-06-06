## FilterOligos.py

Python program to filter hybridization oligos using bwa and primer3. Based on bwa.py and primer3_filter.py from [Chorus2](https://github.com/zhangtaolab/Chorus2) by [zhangtaolab](https://github.com/zhangtaolab).

#### bwa_filter()
Aim: Remove oligos from pipeline with alignment scores (AS) that are too low or suboptimal alignment scores (XS) that are too high.

#### primer_3_filter()
Aim: Remove oligos from pipeline that will make poor hybridization probes because they melt too easily (TM below minimum), are likely to form hairpins (HTM above maximum), or the difference between the two melting temperatures is too small.

### Usage
Program both reads from and outputs to sam format. Uses the [primer3-py](https://pypi.org/project/primer3-py/) library.

Example command:
```
python FilterOligos.py unfiltered.sam filtered.sam
```


## CheckForRepeats.py
Quick-and-dirty script to check for repeats in sam file (dupicate mappings).
