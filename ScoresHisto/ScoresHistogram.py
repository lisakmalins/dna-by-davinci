# 17 June 2019
# Lisa Malins
# ScoresHistogram.py

"""
Counts occurrences of k-mer scores from sam output of CalcKmerScores.py.

Usage:
python ScoresHistogram.py filename.sam

# Specify output file name/path (optional)
python ScoresHistogram.py filename.sam output.txt
"""

from collections import defaultdict
import sys

# Setup file IO
if len(sys.argv) < 2 or len(sys.argv) > 3:
    exit("Usage: python ScoresHistogram.py filename.sam")

try:
    filename = sys.argv[1]
    source = open(filename, 'r')
except:
    exit("File " + str(filename) + " not found")

if len(sys.argv) == 3:
    # Use provided output file name if given
    outputname = sys.argv[2]
else:
    # Otherwise, transform input file name
    outputname = filename.rsplit(".", 1)[0] + "_histo.txt"
output = open(outputname, 'w')

print("Reading scores from", filename)
print("Score histogram will be written to", outputname)


# Get length of file for progress output
source.seek(0,2)
filelength = float(source.tell())
source.seek(0)
percent = 10

# Count scores into default dictionary
scores_dict = defaultdict(int)

line = source.readline()

while(line):
    # Skip header lines
    if line[0] == "@":
        line = source.readline()
        continue
    # Output progress message
    if (source.tell() / filelength * 100) > percent:
        print("Read progress: " + str(percent) + "%")
        percent += 10
    # Get score and put into dictionary
    score = line.split("\t")[-1][5:].rstrip("\n")
    scores_dict[score] += 1

    line = source.readline()

# Output score histogram as CSV
print("Read complete. Writing score histogram to", outputname)

output.write("score, frequency\n")
for s, f in scores_dict.items():
    output.write(str(s) + "," + str(f) + "\n")

print("Write complete. Have a fantastic day!")
