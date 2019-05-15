# 14 May 2019
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

    # Create dummy entry... is this bad?
    counts = {"": {"": {"": 0}}} # Is this bad?

    # Read all k-mers and counts into dictionary
    line = source.readline()
    num_entries = 0
    cur_size = 0;

    while line:
        # Error message for unreadable input
        assert line[0] == ">", \
        "\nUnable to read k-mers and scores due to unexpected input. " + \
        "Line was:\n" + line.rstrip('\n') + "\nfrom " + sourcename

        # Read count and associated sequence
        count = line[1:].rstrip('\n')
        line = source.readline()
        seq = line.rstrip('\n')
        level1 = seq[0:6]
        level2 = seq[0:12]

        # Insert entry

        if not level1 in counts:
            # If there no entry for level1
            counts[level1] = {level2: {seq: count}}
            # print("Making level 1 entry for " + level1)
        elif not level2 in counts[level1]:
            # If there is an entry for level1 but not level2
            counts[level1][level2] = {seq, count}
            # print("Making level 2 entry for " + level2)
        elif not seq in counts[level1][level2]:
            # If there is an entry for level1 and level2
            counts[level1][level2][seq] = count
            # print("Making level 3 entry for " + seq)
        else:
            raise AssertionError("Duplicate entry found for sequence " \
            + seq + " in " + sourcename)

        # # Watch size grow
        # if cur_size != sys.getsizeof(counts):
        #     cur_size = sys.getsizeof(counts)
        #     print("Number of entries = " + str(num_entries) + \
        #     " size in memory = " + str(sys.getsizeof(counts)))

        line = source.readline()
        num_entries += 1

    # Remove dummy entry
    counts.pop("")

    # Output size of dictionary and further commands
    print(str(num_entries) + " kmers and counts read from file " + sourcename)
    print("Current size of dictionary in memory: " + str(sys.getsizeof(counts)))
    print("Further commands: ")
    print("\tcounts # Print all entries in dictionary")
    print("\tnum_entries # Output number of entries in dictionary")
    print("\tsys.getsizeof(counts) # Size of dictionary in bytes")
    # Does not work anymore
    # print("\tcounts['{sequence}'] # Query count for a particular sequence")

except NameError:
    print("No jellyfish dump input file specified.")
    print("Please set filename as follows and try again:\n" + \
    "sourcename = \"{filename.fa}\"")
except FileNotFoundError:
    if sourcename:
        print("File " + sourcename + " not found")
    print("Please set filename as follows and try again:\n" + \
    "sourcename = \"{filename.fa}\"")
