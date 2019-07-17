# 11 July 2019
# Lisa Malins
# CalcKmerScores.py

"""
Reads 17-mers and counts from jellyfish dump file into a nested dictionary.
Then, calculates k-mer scores for 45-mers from either sam file or fasta file.
If oligos are read from fasta, output is also in fasta format
    with the score appended to the header comment.
If oligos are read from sam, output is also in sam format
    with the score appended as KS:i: tag.

Usage:
python CalcKmerScores.py dump.fa oligos.sam scores_output.sam

If you need to calculate scores with 45-mers from multiple files
but using same dictionary, use -i flag to open interactive mode
at close of program and run the following command:
>>> CalcFromSam(nkd, "nextoligofile.sam", "nextoutputfile.sam", log)
"""

import sys
from NestedKmerDict import NestedKmerDict
from time import ctime
try:
    from time import process_time
except:
    from time import clock as process_time #python2
from datetime import timedelta

# not updated, please skip to CalcFromSam
def CalcFromFasta(nkd, oligos, output, dump, log):
    # Setup nested kmer dictionary
    nkd = NestedKmerDict()
    nkd.Populate(dumpname)

    # Setup IO
    source = open(oligoname, 'r')
    output = open(outputname, 'w')

    # Setup log file
    log = open(logname, 'a+')

    # Begin log file with context
    log.write("> Log file for CalcKmerScores.py\n")
    log.write("> Dictionary loaded from jellyfish dump file = " + dumpname + "\n")
    log.write("> 45-mer source file = " + oligoname + "\n")
    log.write("> Output file of 45-mers and k-mer scores = " + outputname + "\n")

    # Read 45-mers and calculate k-mer scores
    header = source.readline().rstrip('\n')
    while header:
        # Error message if line is not a fasta header
        assert header[0] == ">", \
        "\nUnable to read k-mers and scores due to unexpected input. " + \
        "Line was:\n" + line.rstrip('\n') + "\nfrom " + oligofile

        # Read 45-mers and calculate k-mer score
        oligo = source.readline().rstrip('\n')
        # print("Calculating score for " + line) #debug

        # Start with score of zero
        score = 0

        # Loop through 45-mer and query all 17-mers
        for i in range (0, 29):
            try:
                seq = oligo[i:i+17]
                count = nkd.Query(seq)
                score += int(count)

            # If 17-mer not found in dictionary, note in log and skip it
            except:
                log.write("No dictionary entry for " + seq + \
                " from source oligo " + header + "\n")
                continue

        output.write(header + " " + str(score) + "\n")
        output.write(oligo + "\n")
        header = source.readline().rstrip('\n')

    print("Finished writing k-mers and scores in fasta format to " + outputname)
    print("Log written to " + logname)

def CalcFromSam(nkd, oligos, output, log, fast=True):
    if isinstance(oligos, str):
        oligos = open(oligos, 'r')
    if isinstance(output, str):
        output = open(output, 'w')

    time0 = process_time()

    # Begin log file with context
    log = open(log.name, 'a+')
    log.write("> Beginning k-mer score calculation for file = " + oligos.name + " at " + ctime() + "\n")
    log.write("> Output file of 45-mers and k-mer scores = " + output.name + "\n")
    log.write("> Fast mode is " +  "on\n" if fast else "off\n")
    print("Beginning k-mer score calculation for file = " + oligos.name + " at " + ctime())
    print("Output file of 45-mers and k-mer scores = " + output.name)
    print("Fast mode is " +  "on" if fast else "off")

    # Read 45-mers and calculate k-mer scores
    line = oligos.readline()
    while line:
        # Print headers without touching them
        if line[0] == '@':
            output.write(line)
            line = oligos.readline()
            continue

        # Grab oligo and start with score of 0
        oligo = line.split('\t')[9]
        score = 0

        # Loop through 45-mer and query all 17-mers
        for i in range (0, 29):
            try:
                seq = oligo[i:i+17]
                if fast:
                    count = nkd.QueryFast(seq, log)
                else:
                    count = nkd.Query(seq, log)
                score += int(count)

            # If 17-mer not found in dictionary, note in log and skip it
            except:
                log.write("No dictionary entry for " + seq + \
                " from source oligo " + line.split('\t')[0] + "\n")
                continue

        # Write line with k-mer score appended
        output.write(line.rstrip('\n') + "\tKS:i:" + str(score) + "\n")

        line = oligos.readline()


    log.write("> K-mer score calculation for file " + oligos.name + " completed successfully at " + ctime() + "\n")
    log.write("> Scores output at " + output.name + "\n")
    proc_time = process_time() - time0
    log.write("> Calculation time: " + str(timedelta(seconds=proc_time)) + " (total seconds = " + str(proc_time) + ")\n")
    print("Finished writing k-mers and scores in sam format to " + output.name)
    print("Log written to " + log.name)

    oligos.close()
    output.close()


# ----------------main-------------------

usage = "Usage: {dump input file} {oligo input file} {scores output file} " \
"Optional: {custom log file name} {fast mode True/False}"

if len(sys.argv) < 4 or len(sys.argv) > 6:
    exit(usage)

# Verify all files found and not garbage before loading dictionary
# Open jellyfish dump file of 17-mers
try:
    dump = open(sys.argv[1], 'r')
except FileNotFoundError:
    exit("File " + sys.argv[1] + " not found.")
print("Will read counts from " + dump.name)

# Open file of 45-mers
try:
    oligos = open(sys.argv[2], 'r')
except FileNotFoundError:
    exit("File " + sys.argv[2] + " not found.")
print("Will read oligos from " + oligos.name)

# Remember that one time I named the log but forgot to name the output file
# and then it wrote the output and the log in the same place lol that was hilarious
assert sys.argv[3][-3:] != "log", "Make sure you specify an output file\n" + usage
# Open output file
output = open(sys.argv[3], 'w')
print("Will write scores to " + output.name)

# Open log file
logfile = sys.argv[4] if len(sys.argv) > 4 else output.name.split('.', 1)[0] + ".log"
log = open(logfile, 'a+')
print("Logging to " + log.name)
log.write("> Log file for CalcKmerScores.py\n")

# Setup nested kmer dictionary
nkd = NestedKmerDict()
nkd.Populate(dump, log)

if sys.argv[2].split('.')[-1] == "sam":
    CalcFromSam(nkd, oligos, output, log)
elif sys.argv[2].split('.')[-1][:2] == "fa":
    CalcFromFasta(nkd, oligos, output, dump, log)
else:
    exit(usage)
