# 10 July 2019
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
"""

import sys
from NestedKmerDict import NestedKmerDict

# not updated, please skip to CalcFromSam
def CalcFromFasta(oligos, output, dump, log):
    # Setup nested kmer dictionary
    nkd = NestedKmerDict()
    nkd.Populate(dumpname)

    # Setup IO
    source = open(oligoname, 'r')
    output = open(outputname, 'w')

    # Setup log file
    log = open(logname, 'w')

    # Begin log file with context
    log.write("Log file for CalcKmerScores.py\n")
    log.write("Dictionary loaded from jellyfish dump file = " + dumpname + "\n")
    log.write("45-mer source file = " + oligoname + "\n")
    log.write("Output file of 45-mers and k-mer scores = " + outputname + "\n")

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

def CalcFromSam(oligos, output, dump, log):
    # Setup nested kmer dictionary
    nkd = NestedKmerDict()
    nkd.Populate(dump)

    # Begin log file with context
    log.write("Log file for CalcKmerScores.py\n")
    log.write("Dictionary loaded from jellyfish dump file = " + dump.name + "\n")
    log.write("Sam source file = " + oligos.name + "\n")
    log.write("Output file of 45-mers and k-mer scores = " + output.name + "\n")

    # Read 45-mers and calculate k-mer scores
    line = oligos.readline()
    while line:
        # Print headers without touching them
        if line[0] == '@':
            output.write(line)
            line = oligos.readline()
            continue

        # Split line into components and grab oligo
        splitline = line.split('\t')
        oligo = splitline[9]

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
                " from source oligo " + splitline[0] + "\n")
                continue

        # Write line with k-mer score appended
        output.write(line.rstrip('\n') + "\tKS:i:" + str(score) + "\n")

        line = oligos.readline()

    print("Finished writing k-mers and scores in sam format to " + output.name)
    print("Log written to " + log.name)


# ----------------main-------------------

if len(sys.argv) < 3 or len(sys.argv) > 4:
    exit("Usage: {dump input file} {oligo input file} {scores output file}")

# Verify all files found and not garbage before loading dictionary
# Open jellyfish dump file of 17-mers
try:
    dump = open(sys.argv[1], 'r')
except FileNotFoundError:
    exit("File " + sys.argv[1] + " not found.")

# Open file of 45-mers
try:
    oligos = open(sys.argv[2], 'r')
except FileNotFoundError:
    exit("File " + sys.argv[2] + " not found.")

# Open output file
output = open(sys.argv[3], 'w')

# Open log file
try:
    logfile = sys.argv[4]
except:
    logfile = output.name.split('.', 1)[0] + ".log"
log = open(logfile, 'w')


if sys.argv[2].split('.')[-1] == "sam":
    CalcFromSam(oligos, output, dump, log)
elif sys.argv[2].split('.')[-1][:2] == "fa":
    CalcFromFasta(oligos, output, dump, log)
else:
    exit("Usage: {dump input file} {oligo input file} {scores output file}")
