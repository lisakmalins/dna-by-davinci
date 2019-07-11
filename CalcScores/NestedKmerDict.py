# 11 July 2019
# Lisa Malins
# NestedKmerDict.py

"""
Nested k-mer dictionary class which holds 17-mers in 3 levels.
Reads 17-mers from a Jellyfish dump file.
Multiple Jellyfish dump files can be read into same dictionary object.
"""

import sys
from time import ctime
from time import process_time
from datetime import timedelta

class NestedKmerDict():
    def __init__(self):
        self.counts = {"": {"": {"": 0}}}
        self.num_entries = 0
        self.cur_size = 0

    # Read 17-mers from Jellyfish dump file
    # Accepts string of filename or file object
    def Populate(self, source, log=open("/dev/fd/1", 'w')):
        time0 = process_time()

        # If string of filename passed, reassign variable to be file object
        if isinstance(source, str):
            try:
                source = open(source, 'r')
            except FileNotFoundError:
                exit("File " + source.name + " not found")

        # Read all k-mers and counts into dictionary
        print("Reading kmer counts from file " + source.name + "...")
        log.write("Kmer loading from " + source.name + " began at time " + ctime() + "\n")
        line = source.readline()

        while line:
            # Error message for unreadable input
            assert line[0] == ">", \
            "\nUnable to read k-mers and scores due to unexpected input. " + \
            "Line was:\n" + line.rstrip('\n') + "\nfrom " + source.name

            # Read count and associated sequence
            count = line[1:].rstrip('\n')

            line = source.readline()
            seq = line.rstrip('\n')

            level1 = seq[0:6]
            level2 = seq[6:12]
            level3 = seq[12:17]

            # Insert entry

            # If there no entry for level1
            if not level1 in self.counts:
                self.counts[level1] = {level2: {level3: count}}

            # If there is an entry for level1 but not level2
            elif not level2 in self.counts[level1]:
                self.counts[level1][level2] = {level3: count}

            # If there is an entry for level1 and level2
            elif not level3 in self.counts[level1][level2]:
                self.counts[level1][level2][level3] = count

            # If entry already exists, raise error (should be no duplicates in file)
            else:
                raise AssertionError("Duplicate entry found for sequence " \
                + seq + " in " + source.name)

            self.num_entries += 1

            # Watch size grow
            # if cur_size != sizeofdict(counts):
            #     cur_size = sizeofdict(counts)
            #     print("Number of entries =", num_entries, \
            #     "\ttop level size =", sys.getsizeof(counts), \
            #     "\ttotal size =", cur_size)

            line = source.readline()

        # Close source file and remove dummy entry if necessary
        source.close()
        self.counts.pop("") if "" in self.counts else None

        # Output size of dictionary
        proc_time = process_time() - time0
        print(str(self.num_entries) + " kmers and counts read from file " + source.name)
        log.write("Kmer loading from " + source.name + " completed at time " + ctime() + "\n")
        log.write("Load time: " + str(timedelta(seconds=proc_time)) + " (total seconds = " + str(proc_time) + ")\n")
        log.write("Total size in memory is " + str(self.Size()) + "\n")

        # Print help if running in interactive mode
        if not sys.argv[0]:
            print("\nFurther commands:")
            self.Help()

    # Converts sequence to reverse complement
    def RC(self, seq):
        rc = ""
        for letter in seq:
            if letter == 'A':
                rc += 'T'
            elif letter == 'C':
                rc += 'G'
            elif letter == 'G':
                rc += 'C'
            elif letter == 'T':
                rc += 'A'
            else:
                print("LOL that's literally not in my vocabulary")
        return rc

    # Find count for k-mer or its reverse complement
    # Always checks both forward and reverse and logs if both are found
    def Query(self, seq, log=open("/dev/fd/1", 'w')):
        matches = 0
        try:
            fcount = self.counts[seq[0:6]][seq[6:12]][seq[12:17]]
            matches += 1
        except KeyError:
            fcount = 0
        try:
            rc = self.RC(seq)
            rcount = self.counts[rc[0:6]][rc[6:12]][rc[12:17]]
            matches += 1
        except KeyError:
            rcount = 0
        except TypeError:
            print("Please enter a DNA sequence in quotes")

        if matches == 1:
            return fcount or rcount
        elif matches == 2:
            log.write("Both " + seq + " and reverse complement " + rc + " found in dictionary\n")
            return max(fcount, rcount)
        else:
            log.write(seq + " not found in dictionary\n")
            return 0

    # Find count for k-mer or its reverse complement
    # Only checks for reverse complement if forward not found
    def QueryFast(self, seq, log=open("/dev/fd/1", 'w')):
        try:
            return self.counts[seq[0:6]][seq[6:12]][seq[12:17]]
        except KeyError:
            rc = self.RC(seq)
            try:
                return self.counts[rc[0:6]][rc[6:12]][rc[12:17]]
            except KeyError:
                log.write(seq + " not found in dictionary\n")

    def Size(self, sum=0, verbose=False):
        return self._Size(self.counts, sum, verbose)

    def _Size(self, d, sum=0, verbose=False):
        # Execute this block for innermost values only
        # Add size of final value and stop recursion
        if not isinstance(d, dict):
            # Add size of final value
            sum += sys.getsizeof(d)
            return sum

        # Execute this block for dictionaries
        # For each item, add size of key and recur on members
        for item in d:
            # Add size of string key
            sum += sys.getsizeof(item)
            # Recur on members of dict
            sum = self._Size(d[item], sum)

        # Add size of dictionary
        sum += sys.getsizeof(d)
        return sum

    def PrintAll(self):
        return self.counts
    def NumEntries(self):
        return self.num_entries

    # Help function designed for interactive mode
    def Help(self):
        print("# Print all entries in dictionary")
        print("PrintAll()")
        print("# Output number of entries in dictionary")
        print("NumEntries()")
        print("# Output size of dictionary in bytes")
        print("Size()")
        print("# Query count for a particular sequence")
        print("Query(sequence)")
        print("# Display this help menu")
        print("Help()")

# Demo
if __name__ == '__main__':
    nkd = NestedKmerDict()
    nkd.Populate("fakedump.fa")
    # nkd.Populate("monkeywrench.fa")
    # nkd.Query("GGGGGGGGGGGGGGGGG")
    print("Size =", nkd.Size())
