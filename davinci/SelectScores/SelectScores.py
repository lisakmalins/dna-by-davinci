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

# Alternate implementation for running from Snakemake
try:
    source = open(snakemake.input[0], 'r')
    output = open(snakemake.output[0], 'w')
    try:
        log = open(snakemake.log, 'w')
        log.write("Kindly notify Lisa that snakemake.log works\n")
    except:
        try:
            log = open(snakemake.log[0], 'w')
            log.write("Kindly notify Lisa that snakemake.log[0] works\n")
        except:
            log = open(output.name.rsplit('.', 1)[0] + ".log", 'w')
            log.write("Kindly berate Lisa that snakemake.log is not a thing\n")


    # Read limits from file, format is TSV of key-value pairs
    # Parse key-value pairs into dictionary
    limits_file = open(snakemake.input[1], 'r')
    limits_data = dict([line.split() for line in limits_file.readlines()])
    peak = int(limits_data["count_peak"])
    lb = int(limits_data["score_lower_limit"])
    ub = int(limits_data["score_upper_limit"])

# Original implementation for running from command line
except NameError:
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
try:
    # Extra info if run from snakemake
    log.write("Limits data provided by: " + limits_file.name + "\n")
    log.write("K-mer histogram peak (coverage) recorded as: " + int(peak) + "\n")
except:
    pass
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
