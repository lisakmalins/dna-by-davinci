# 14 May 2019
# Lisa Malins
# LoadKmerDict.py

"""
Reads k-mers and counts from Jellyfish dump file into a dictionary
in order to measure how much memory is required to hold them.

From python shell:
dumpfile = "{filename.fa}"
exec(open("LoadKmerDict.py").read())
"""

import sys

try:
    # Make sure script is being run in interactive mode
    assert sys.argv[0] != "LoadKmerDict.py", "\n" + \
    "You ran LoadKmerDict.py in script mode.\n" + \
    "Please open the python shell and rerun in interactive mode " + \
    "using the following commands:\n" + \
    "dumpfile = \"{filename.fa}\"\n" + \
    "exec(open(\"LoadKmerDict.py\").read())"

    # Open file
    source = open(dumpfile, 'r')
    print("Reading kmer counts from file " + dumpfile + "...")

    # Create dummy entry
    counts = {"": {"": {"": 0}}}

    # Read all k-mers and counts into dictionary
    line = source.readline()
    num_entries = 0
    cur_size = 0;

    while line:
        # Error message for unreadable input
        assert line[0] == ">", \
        "\nUnable to read k-mers and scores due to unexpected input. " + \
        "Line was:\n" + line.rstrip('\n') + "\nfrom " + dumpfile

        # Read count and associated sequence
        count = line[1:].rstrip('\n')

        line = source.readline()
        seq = line.rstrip('\n')

        level1 = seq[0:6]
        level2 = seq[0:12]

        # Insert entry

        # If there no entry for level1
        if not level1 in counts:
            counts[level1] = {level2: {seq: count}}

        # If there is an entry for level1 but not level2
        elif not level2 in counts[level1]:
            counts[level1][level2] = {seq, count}

        # If there is an entry for level1 and level2
        elif not seq in counts[level1][level2]:
            counts[level1][level2][seq] = count

        # If entry already exists, raise error (should be no duplicates in file)
        else:
            raise AssertionError("Duplicate entry found for sequence " \
            + seq + " in " + dumpfile)

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
    print(str(num_entries) + " kmers and counts read from file " + dumpfile)
    print("Current size of dictionary in memory: " + str(sys.getsizeof(counts)))
    print("Further commands: ")
    print("\t# Print all entries in dictionary")
    print("\tcounts")
    print("\t# Output number of entries in dictionary")
    print("\tnum_entries")
    print("\t# Output size of dictionary in bytes")
    print("\tsys.getsizeof(counts)")
    print("\t# Query count for a particular sequence")
    print("\tcounts['{first 6 letters}']['{first 12 letters}']" + \
    "['{all 17 letters}']")

except NameError:
    print("No jellyfish dump input file specified.")
    print("Please set filename as follows and try again:\n" + \
    "dumpfile = \"{filename.fa}\"")
except FileNotFoundError:
    if dumpfile:
        print("File " + dumpfile + " not found")
    print("Please set filename as follows and try again:\n" + \
    "dumpfile = \"{filename.fa}\"")
