# 11 July 2019
# Lisa Malins
# NestedKmerDict.py

"""
Nested k-mer dictionary class which holds 17-mers in 3 levels.
Reads 17-mers from a Jellyfish dump file.
Multiple Jellyfish dump files can be read into same dictionary object.
"""

import sys
import gc
from time import ctime
try:
    from time import process_time
except:
    from time import clock as process_time #python2
from datetime import timedelta

class NestedKmerDict():
    def __init__(self, source=None):
        # Empty counts dictionary
        self.counts = {"": {"": {"": 0}}}

        # Size info
        self.num_entries = 0
        self.cur_size = 0

        # Whether duplicates found (for logging)
        self.dup_found = False

        # Percent cleared (for status messages in clear)
        self.percent_emptied = 0

        # Deluxe constructor: Can populate at same time
        if source is not None:
            self.Populate(source)

    def __del__(self):
        sys.stderr.write("\nNkd destructor says: Hey I'm running\n")

        self.Clear()

        try:
            sys.stderr.write("Nkd destructor says: Current size: {}\n".format(str(self.Size())))
            # sys.stderr.write("Nkd destructor says: Current contents: \n")
            # sys.stderr.write(str(self.PrintAll()) + "\n")
        except:
            sys.stderr.write("Nkd destructor says: Could not display current size (EXPECTED)\n")

        try:
            gc.collect()
            sys.stderr.write("Nkd destructor says: Garbage collected\n")
        except:
            sys.stderr.write("Nkd destructor says: Could not garbage collect\n")

        try:
            del self.counts
            sys.stderr.write("Nkd destructor says: Deleted self.counts\n")
        except:
            sys.stderr.write("Nkd destructor says: Could not delete self.counts\n")

        try:
            gc.collect()
            sys.stderr.write("Nkd destructor says: Garbage collected\n")
        except:
            sys.stderr.write("Nkd destructor says: Could not garbage collect\n")
        sys.stderr.write("Nkd destructor says: Nested kmer dict dying\n")

    # Recursively empty dictionary and print helpful status messages
    def Clear(self, verbose=False):
        sys.stderr.write("\nNkd clear says: Hey what's up I'm your friendly neighborhood clear function\n")
        sys.stderr.write("Nkd clear says: Let me calculate how much garbage I have to collect\n")
        sys.stderr.write("Nkd clear says: I have {} entries slash {} bytes of memory to clear so feel free to grab a drink...\n".format(self.num_entries, str(self.Size())))

        self._Clear(self.counts, self.num_entries)

        sys.stderr.write("Nkd clear says: My work here is done.\t{}\n".format(ctime()))


    def _Clear(self, obj, orig_num_entries, level=1, verbose=False, hella_verbose=False, debug_percent=False):
        # Status counter
        if debug_percent:
            print("\nNum entries: ", self.num_entries)
            print("percent complete: ", str((1 - self.num_entries / orig_num_entries) * 100))
            print("Current % counter: " + str(self.percent_emptied))
        if ((1 - self.num_entries / orig_num_entries) * 100 > self.percent_emptied):
            # May need to increment more than once for extremely small datasets
            while ((1 - self.num_entries / orig_num_entries) * 100 > self.percent_emptied + 9.99):
                self.percent_emptied += 10

            sys.stderr.write("Garbage collection {}% complete\t\t {}\n".format(\
            str(self.percent_emptied), ctime()))
            self.percent_emptied += 10

        # When terminal key found, decrement num_entries and return
        if not isinstance(obj, dict):
            if verbose:
                sys.stderr.write("Terminal item found:" + str(obj) + "\n")
            self.num_entries -= 1
            return

        # If current object is a dictionary, recur on all items
        # After all nested items have returned, clear this level and return
        if verbose:
            print("Nested dict found:")
            if hella_verbose:
                print(obj)

        for item in obj:
            if verbose:
                sys.stderr.write("Recursing on {}\n".format(item))
            self._Clear(obj[item], orig_num_entries, level + 1, verbose, hella_verbose)

        # Clear + delete all inner objects
        if level > 1:
            obj.clear()
            assert len(obj) == 0, "Length of object is not 0 after clearing"
            del obj
            gc.collect()
        # Extra debug for top-level
        else:
            sys.stderr.write("FINAL BOSS: top-level reached\n")
            num_top_items = len(obj)
            top_progress = 10
            sys.stderr.write("{} top-level items to remove\n".format(num_top_items))

            # Pop all items
            while len(obj):
                obj.popitem()

                # Item pop progress counter
                if ( (1 - (len(obj) / num_top_items)) * 100 > top_progress):
                    sys.stderr.write("Top-level item removal {} percent complete\t{}\n".format(\
                    str(top_progress), ctime()))
                    top_progress += 10
                    gc.collect()

            assert len(obj) == 0, "Length of top level is not 0 after clearing"
            sys.stderr.write("Done popping items at {}\n".format(ctime()))
            del obj

        # Extra stats when num_entries nears zero
        if self.num_entries < 5 and level <= 2:
            sys.stderr.write(" ".join(["level", str(level), "top-level-length/size/num_entries:\t"]))
            sys.stderr.write(" / ".join([ str(len(self.counts)), str(self.Size()), str(self.num_entries) + "\n" ] ) )

        # Verify empty before return
        if level == 1:
            sys.stderr.write("Immabout to return, here are my contents and top-level length:\n")
            print(self.PrintAll())
            print(len(self.counts))


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
        sys.stderr.write("Reading kmer counts from file " + source.name + "...\n")
        sys.stderr.write("Logging to " + log.name + "\n")
        log.write("Kmer loading from " + source.name + " began at time " + ctime() + "\n")
        log.flush()
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
        sys.stderr.write(str(self.num_entries) + " kmers and counts read from file " + source.name + "\n")
        sys.stderr.write("Calculating memory size...\n")
        self.cur_size = self.Size()
        sys.stderr.write("Memory size is " + str(self.cur_size) + " bytes.\n")
        log.write("Kmer loading from " + source.name + " completed at time " + ctime() + "\n")
        log.write("Load time: " + str(timedelta(seconds=proc_time)) + " (total seconds = " + str(proc_time) + ")\n")
        log.write("Total size in memory is " + str(self.cur_size) + " bytes for " + str(self.num_entries) + " entries\n")
        log.flush()

        # Print help if running in interactive mode
        if not sys.argv[0]:
            print("\nFurther commands:")
            self.Help()

    # Converts sequence to reverse complement
    def RC(self, seq):
        rc = ""
        # syntax [::-1] reverses string
        for letter in seq[::-1]:
            if letter == 'A':
                rc += 'T'
            elif letter == 'C':
                rc += 'G'
            elif letter == 'G':
                rc += 'C'
            elif letter == 'T':
                rc += 'A'
            else:
                raise ValueError("Function RC only accepts A, C, G, and T")
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
            sys.stderr.write("Please enter a DNA sequence in quotes\n")

        if matches == 1:
            return fcount or rcount
        elif matches == 2:
            log.write("Both " + seq + " and reverse complement " + rc + " found in dictionary\n")
            return max(fcount, rcount)
        else:
            raise KeyError

    # Find count for k-mer or its reverse complement
    # Only checks for reverse complement if forward not found
    def QueryFast(self, seq, log=open("/dev/fd/1", 'w')):
        try:
            return self.counts[seq[0:6]][seq[6:12]][seq[12:17]]
        except KeyError:
            rc = self.RC(seq)
            return self.counts[rc[0:6]][rc[6:12]][rc[12:17]]


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
