# 11 July 2019
# Lisa Malins
# CalcKmerScores.py

"""
Reads 17-mers and counts from jellyfish dump file into a nested dictionary.
Then, calculates k-mer scores for 45-mers in sam file.
Output is also in sam format with the score appended as KS:i: tag.

Usage:
python CalcKmerScores.py dump.fa oligos.sam scores_output.sam
    Optional: custom.log {fast mode True/False}

If you need to calculate scores with 45-mers from multiple files
but using same dictionary, use -i flag to open interactive mode
at close of program and run the following command:
>>> CalcFromSam(nkd, "nextoligofile.sam", "nextoutputfile.sam", log)
"""

import sys
import gc
from NestedKmerDict import NestedKmerDict
from time import ctime
try:
    from time import process_time
except:
    from time import clock as process_time #python2
from datetime import timedelta


# Calculate k-mer scores of oligos from sam file.
# Needs nested k-mer dictionary object, oligo source file, and output file
# fast=False will check both forward and reverse k-mers and log if both are found.
# fast=True will only check for reverse if forward not found.
# log_missing should be either False or an output file object.
def CalcFromSam(nkd, oligos, output, log, fast=True, log_missing=False):
    if isinstance(oligos, str):
        oligos = open(oligos, 'r')
    if isinstance(output, str):
        output = open(output, 'w')

    time0 = process_time()

    # Begin log file with context
    log = open(log.name, 'a')
    log.write("Beginning k-mer score calculation for oligo file " + oligos.name + " at " + ctime() + "\n")
    log.write("Output file of 45-mers and k-mer scores = " + output.name + "\n")
    log.write("Fast mode is " + ("on\n" if fast else "off\n"))
    print("Beginning k-mer score calculation for file = " + oligos.name + " at " + ctime())
    print("Output file of 45-mers and k-mer scores = " + output.name)
    print("Fast mode is " + ("on" if fast else "off"))

    # Read 45-mers and calculate k-mer scores
    num_missing = 0
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
                num_missing += 1
                if log_missing:
                    log_missing.write("No dictionary entry for " + seq + \
                    " from source oligo " + line.split('\t')[0] + "\n")
                continue

        # Write line with k-mer score appended
        output.write(line.rstrip('\n') + "\tKS:i:" + str(score) + "\n")

        line = oligos.readline()

    proc_time = process_time() - time0
    msg = "Calculation time: " + str(timedelta(seconds=proc_time)) + " (total seconds = " + str(proc_time) + ")\n"
    log.write("K-mer score calculation for file " + oligos.name + " completed successfully at " + ctime() + "\n")
    log.write(msg)
    log.write("Scores output at " + output.name + "\n")
    log.write("{} k-mers not found in dictionary\n".format(str(num_missing)))
    sys.stderr.write("Finished writing k-mers and scores in sam format to " + output.name + " at " + ctime() + "\n")
    sys.stderr.write(msg)
    sys.stderr.write("Log written to " + log.name + "\n")
    sys.stderr.write("{} k-mers not found in dictionary\n".format(str(num_missing)))
    if log_missing:
        log.write("Missing k-mers written to " + missing.name + "\n")
        sys.stderr.write("Missing k-mers written to " + missing.name + "\n")

    del nkd
    gc.collect()
    return

# ----------------main-------------------

if __name__ == "__main__":

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

    # Open main log file
    logfile = sys.argv[4] if len(sys.argv) > 4 else output.name.rsplit('.', 1)[0] + ".log"
    log = open(logfile, 'w')
    print("Logging to " + log.name)
    log.write("Log file for CalcKmerScores.py\n")
    log.flush()

    # Separate log file for missing k-mers
    log_missing = True
    if log_missing:
        missing = open(log.name + ".missing", 'w')


    # Setup nested kmer dictionary
    nkd = NestedKmerDict()
    nkd.Populate(dump, log)
    dump.close()

    CalcFromSam(nkd, oligos, output, log, fast=True, log_missing=missing)
