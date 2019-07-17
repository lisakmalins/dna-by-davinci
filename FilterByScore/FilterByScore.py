# 16 July 2019
# Lisa Malins
# FilterByScore.py

"""
Filters output of CalcKmerScores.py by k-mer score.

Usage (output filename is optional):
python FilterByScore.py {input filename} {lower bound} {upper bound} {output filename}
"""

import sys

# Setup file IO
usage = "Usage: python FilterByScore.py {input filename} {lower bound} {upper bound}"

if len(sys.argv) < 4 or len(sys.argv) > 5:
    exit(usage)

try:
    source = open(sys.argv[1], 'r')
except FileNotFoundError:
    if not sys.argv[1].isalpha():
        exit(usage)
    else:
        exit("File " + str(sys.argv[1]) + " not found")

try:
    output = open(sys.argv[4], 'w')
except IndexError:
    output = open(sys.argv[1].rsplit(".", 1)[0] + "_KS_filtered.sam", 'w')

if not sys.argv[2].isdigit() and sys.argv[3].isdigit():
    exit(usage)

lb, ub = int(sys.argv[2]), int(sys.argv[3])

line = source.readline()

while line:
    # Print all headers
    if line[0] == "@":
        output.write(line)
        line = source.readline()
        continue

    # Get k-mer score
    score = int(line.split('\t')[15].rstrip('\n')[5:])

    # Output lines with k-mer scores in range
    if score in range (lb, ub):
        output.write(line)

    # # Debug: print rejected lines
    # else:
    #     print(line.rstrip('\n'), "rejected")

    line = source.readline()
