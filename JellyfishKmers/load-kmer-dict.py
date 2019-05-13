# 13 May 2019
# Lisa Malins
# load-kmer-dict.py

"""
Reads k-mers and counts from Jellyfish dump file into a dictionary
in order to measure how much memory is required to hold them.

From python shell:
sourcename = "{filename.fa}"
exec(open("load-kmer-dict.py").read())
"""

import sys

try:
    # Make sure script is being run in interactive mode
    assert sys.argv[0] != "load-kmer-dict.py", "\n" + \
    "You ran load-kmer-dict.py in script mode.\n" + \
    "Please open the python shell and rerun in interactive mode " + \
    "using the following commands:\n" + \
    "sourcename = \"{filename.fa}\"\n" + \
    "exec(open(\"load-kmer-dict.py\").read())"

    # Open file
    source = open(sourcename, 'r')
    print("Reading kmer counts from file " + sourcename + "...")

    # Read all k-mers and counts into dictionary
    line = source.readline()
    counts = dict()
    while line:
        assert line[0] == ">", \
        "\nUnable to read k-mers and scores due to unexpected input. " + \
        "Line was:\n" + line.rstrip('\n') + "\nfrom " + sourcename
        count = line[1:].rstrip('\n')
        line = source.readline()
        seq = line.rstrip('\n')
        counts[seq] = count
        line = source.readline()

    # Output size of dictionary and further commands
    print(str(len(counts)) + " kmers and counts read from file " + sourcename)
    print("Current size of dictionary in memory: " + str(sys.getsizeof(counts)))
    print("Further commands: ")
    print("\tEnter counts to print all entries in dictionary")
    print("\tEnter sys.getsizeof(counts) to find size of dictionary")
    print("\tEnter counts['sequence'] to see count for a particular sequence")

except NameError:
    print("No jellyfish dump input file specified.")
    print("Please set filename as follows and try again:\n" + \
    "sourcename = \"{filename.fa}\"")
except FileNotFoundError:
    if sourcename:
        print("File " + sourcename + " not found")
    print("Please set filename as follows and try again:\n" + \
    "sourcename = \"{filename.fa}\"")
