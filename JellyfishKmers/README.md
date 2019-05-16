LoadKmerDict.py reads k-mers and counts from Jellyfish dump file into a dictionary in order to measure how much memory is required to hold them.

dump100.fa is a tiny Jellyfish dump file for testing. It is the first 100 lines of a real Jellyfish dump file for maize.

From python shell:

`dumpfile = "{filename.fa}"`

`exec(open("LoadKmerDict.py").read())`
