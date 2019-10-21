# 16 July 2019
# Lisa Malins
# SelectScores.py

"""
Filters output of CalcKmerScores.py by k-mer score.

Usage (output and log filenames are optional):
python SelectScores.py {input filename} {lower bound} {upper bound} {output filename} {log filename}
"""

import sys
from time import ctime
try:
    from time import process_time
except ImportError:
    from time import clock as process_time #python2
from datetime import timedelta

# Setup file IO
usage = "Usage: python SelectScores.py {input filename} {lower bound} {upper bound}"

# Verify number of arguments
if len(sys.argv) < 4 or len(sys.argv) > 6:
    exit(usage)

# Verify input filename
try:
    source = open(sys.argv[1], 'r')
except FileNotFoundError:
    if not sys.argv[1].isalpha():
        exit(usage)
    else:
        exit("File " + str(sys.argv[1]) + " not found")

# Use default if no output filename provided
try:
    output = open(sys.argv[4], 'w')
except IndexError:
    output = open(sys.argv[1].rsplit('.', 1)[0] + "_KS_filtered.sam", 'w')

# Verify upper and lower bounds
if not sys.argv[2].lstrip('-').isdigit() and sys.argv[3].lstrip('-').isdigit():
    print("alice")
    exit(usage)
try:
    lb, ub = int(sys.argv[2]), int(sys.argv[3])
    if lb < 0 or ub < 0: raise ValueError
except ValueError:
    exit("Please provide score ranges as positive integers\n" + usage)

# Setup log file
try:
    log = open(sys.argv[5], 'w')
except IndexError:
    log = open(output.name.rsplit('.', 1)[0] + ".log", 'w')

# Echo arguments to screen
print("Will read oligos from " + source.name)
print("Will write scores in range (" + str(lb) + ", " + str(ub) + ") to " + output.name)
print("Will log to " + log.name)

# Setup file progress counter
print("Scanning file length...")
source.seek(0,2)
filelength = float(source.tell())
source.seek(0)
percent = 10

# Begin log with context
log.write("Log file for SelectScores.py\n")
log.write("Oligos and scores read from: " + source.name + "\n")
log.write("File size: " + str(filelength) + " bytes\n")
log.write("Oligos with scores in range (" + str(lb) + ", " + str(ub) + ") wrtten to: " + output.name + "\n")
log.write("Filtering began: " + ctime() + "\n")

print("Filtering oligos...")

while True:
    line = source.readline()
    if not line: break

    # Print progress to screen
    if (source.tell() / filelength * 100) > percent:
        print("Progress: " + str(percent) + "%")
        percent += 10

    # Write all headers
    if line[0] == "@":
        output.write(line)
        continue

    # Get k-mer score
    score = int(line.split('\t')[15].rstrip('\n')[5:])

    # Output lines with k-mer scores in range
    if score in range (lb, ub):
        output.write(line)

    # # Debug: print rejected lines
    # else:
    #     print(line.rstrip('\n'), "rejected")

log.write("Filtering completed successfully: " + ctime() + "\n")
log.write("Total time: " + str(timedelta(seconds=process_time())) + "\n")
print("Write completed successfully. Filtered scores written to " + output.name)
print("Log written to " + log.name)
