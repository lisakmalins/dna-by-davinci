## ScoresHistogram.py

Counts occurrences of k-mer scores from sam output of CalcKmerScores.py.

### Usage:
```
python ScoresHistogram.py filename.sam

# Default output name "{filename}_histo.txt" in same directory

# Specify output file name/path (optional)
python ScoresHistogram.py filename.sam outputname.txt
```

Output is a text file beginning with `score, frequency` followed by a histogram of comma-separated score-frequency pairs. The file can be read into R to make a visual histogram.
