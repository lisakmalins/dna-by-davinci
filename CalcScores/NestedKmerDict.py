# 14 May 2019
# Lisa Malins
# LoadKmerDict.py

"""
Reads k-mers and counts from Jellyfish dump file into a dictionary
in order to measure how much memory is required to hold them.

From python shell:
dumpfile = "{filename.fa}"
exec(open("LoadKmerDict.py").read())

from NestedKmerDict import NestedKmerDict
nkd = NestedKmerDict()
"""

import sys
from GetSizeOfDict import sizeofdict

class NestedKmerDict():
    def __init__(self):
        self.counts = {"": {"": {"": 0}}}
        self.num_entries = 0
        self.cur_size = 0

    def Populate(self, dumpfile):
        try:
            # Open file
            source = open(dumpfile, 'r')
            print("Reading kmer counts from file " + dumpfile + "...")

            # Read all k-mers and counts into dictionary
            line = source.readline()

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
                if not level1 in self.counts:
                    self.counts[level1] = {level2: {seq: count}}

                # If there is an entry for level1 but not level2
                elif not level2 in self.counts[level1]:
                    self.counts[level1][level2] = {seq: count}

                # If there is an entry for level1 and level2
                elif not seq in self.counts[level1][level2]:
                    self.counts[level1][level2][seq] = count

                # If entry already exists, raise error (should be no duplicates in file)
                else:
                    raise AssertionError("Duplicate entry found for sequence " \
                    + seq + " in " + dumpfile)

                self.num_entries += 1

                # Watch size grow
                # if cur_size != sizeofdict(counts):
                #     cur_size = sizeofdict(counts)
                #     print("Number of entries =", num_entries, \
                #     "\ttop level size =", sys.getsizeof(counts), \
                #     "\ttotal size =", cur_size)

                line = source.readline()

            # Close source file and remove dummy entry
            source.close()
            self.counts.pop("")

            # Output size of dictionary
            print(str(self.num_entries) + " kmers and counts read from file " + dumpfile)

            # Print help if running in interactive mode
            if not sys.argv[0]:
                print("\nFurther commands:")
                self.Help()

        except FileNotFoundError:
            if dumpfile:
                print("File " + dumpfile + " not found")
            print("Please set filename as follows and try again:\n" + \
            "dumpfile = \"{filename.fa}\"")


    def PrintAll(self):
        return self.counts
    def NumEntries(self):
        return self.num_entries
    def Size(self):
        cur_size = sizeofdict(self.counts)
        return cur_size

    # Help function designed for interactive mode
    def Help(self):
        print("# Print all entries in dictionary")
        print("PrintAll()")
        print("# Output number of entries in dictionary")
        print("NumEntries()")
        print("# Output size of dictionary in bytes")
        print("Size()")
        print("# Query count for a particular sequence")
        print("counts['{first 6 letters}']['{first 12 letters}']" + \
        "['{all 17 letters}']")
        print("# Display this help menu")
        print("Help()")
