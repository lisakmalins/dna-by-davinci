# 29 July 2019
# Lisa Malins
# GetOligos.py

"""
Script slices fasta genome assembly into overlapping oligos.

Required argument: genome filename
Optional arguments: kmer size, step size, output filename, log filename, sequences to keep.

For more usage information:
python GetOligos.py --help
"""

import sys
import argparse
from os import stat
from time import ctime
try:
    from time import process_time
except:
    from time import clock as process_time #python2
from datetime import timedelta

class HeaderException(Exception):
    pass

class Kmer:

    # Starts a new k-mer with given index and sequence
    def StartNew(self, id, seq):
        self.seq = seq
        self.id = id
        self.index = 1

    # Appends next sequence and forgets the same # characters that is added
    def Advance(self, addition):
        self.seq = self.seq[len(addition):] + addition
        self.index += len(addition)


# Finds headers returns the id.
# If user specified sequences to read, this function skips over sequences not in user's list.
def NextHeader(source, log, line=" ", read_all_seqs=True, seqs_to_read=""):
    while True:
        # If the next line is a header
        if line[0] == ">":
            id = line[1:].split()[0]
            if read_all_seqs:
                return id
            elif id in seqs_to_read:
                return id
            else:
                log.write("Ignored sequence:\n" + line)

        line = args.genome.readline()

# Reads the specified number of characters from source and returns string
def ReadChars(source, length):
    addition = ""
    while len(addition) < length:
        next_letter = source.read(1)

        # Append letters
        if next_letter.isalpha():
            addition += next_letter
        # Ignore whitespace
        elif next_letter.isspace():
            continue
        # Throw headers back
        elif next_letter == '>':
            line = ">" + source.readline()
            raise HeaderException(line)
        # Throw EOF when out of letters
        elif not next_letter:
            raise EOFError
        # Anything else is a problem
        else:
            print("Unexpected character", next_letter, "from file", source.name)
            sys.exit(1)
    return addition


def read_args():
    parser = argparse.ArgumentParser("Get oligos from a genome assembly.\n")
    parser.add_argument("-g", "--genome", type=argparse.FileType('r'), required=True, help="filename of genome assembly to get oligos from")
    parser.add_argument("-m", "--mer-size", type=int, default=45, help="desired oligo size in bases")
    parser.add_argument("-s", "--step-size", type=int, default=3, help="number of bases between start of consecutive oligos")
    parser.add_argument("-o", "--output", type=argparse.FileType('w'), help="output filename")
    parser.add_argument("-l", "--log", type=argparse.FileType('w'), help="log filename")

    sequence_args = parser.add_mutually_exclusive_group()
    sequence_args.add_argument("--sequences", type=str, nargs="+", help="space-separated list of sequences to get oligos from (default: all)")
    sequence_args.add_argument("--seqfile", type=argparse.FileType('r'), help="file with list of sequences to get oligos from, one per line (default: all)")

    args = parser.parse_args()

    # Validate mer size and step size
    if args.mer_size == 0:
        raise Exception("Error: mer size must be greater than 0")
    if args.step_size == 0:
        raise Exception("Error: step size must be greater than 0")

    # Set default output and log filenames
    if args.output is None:
        args.output = open(args.genome.name.rsplit('.', 1)[0] + "_" + str(args.mer_size) + "mers.fa", 'w')
    if args.log is None:
        args.log = open(args.output.name.rsplit('.', 1)[0] + ".log", 'w')

    return args

#-------------------main-----------------------
args = read_args()

# Read all sequences unless user specified which to read
read_all_seqs = True
if args.sequences:
    read_all_seqs = False
    seqs_to_read = args.sequences
elif args.seqfile:
    read_all_seqs = False
    seqs_to_read = list(map(lambda x: str(x).strip(), args.seqfile.readlines()))
    args.seqfile.close()


# Begin log file
args.log.write("Log file for GetOligos.py\n")
args.log.write("Genome file: " + args.genome.name + "\n")
if not read_all_seqs:
    args.log.write("Reading the following sequences only:\n")
    args.log.write("\n".join(seqs_to_read))
    args.log.write("\n")
args.log.write("Oligo size: " + str(args.mer_size) + "\t")
args.log.write(" Step size: " + str(args.step_size) + "\n")
args.log.write("Oligos written to: " + args.output.name + "\n\n")


# Get length of file for progress output
filelength = float(stat(args.genome.name).st_size)
percent = 10

# Begin status messages to screen
print("Reading " + str(args.mer_size) + "-mers with step size of " + str(args.step_size) + \
" from " + args.genome.name + " and writing to " + args.output.name)
print("Logging to: " + args.log.name)
print("Genome slicing into oligos beginning at " + ctime())
print("Progress messages will be displayed here.")
print("Progress messages may not reach 100% because of ignored sequences.")
args.log.write("\nGenome slicing into oligos beginning at " + ctime() + "\n\n")

time0 = process_time()

# Create Kmer object
kmer = Kmer()
line = " "

# Main loop for getting oligos
try:
    while True:
        # Progress messages
        if (args.genome.tell() / filelength * 100 >= percent):
            print("Read progress : " + str(percent) + "%")
            percent += 10

        # Get header, start new k-mer
        if read_all_seqs:
            id = NextHeader(args.genome, args.log, line)
        else:
            id = NextHeader(args.genome, args.log, line, read_all_seqs=False, seqs_to_read=seqs_to_read)

        kmer.StartNew(id, ReadChars(args.genome, 45))

        # Get k-mers until header encountered
        while True:
            args.output.write(">" + str(id) + "_" + str(kmer.index) + "\n")
            args.output.write(kmer.seq + "\n")
            try:
                kmer.Advance(ReadChars(args.genome, 3))
            # Catch header line and carry to next main loop iteration
            except HeaderException as e:
                line = e.args[0]
                break

# Aaaaaaand stick the landing
except (IndexError, EOFError) as e:
    args.genome.close()

    proc_time = process_time() - time0

    print("Finished writing " + str(args.mer_size) + "-mers to " + args.output.name)
    print("Program finished successfully at " + ctime())
    print("Total time " + str(timedelta(seconds=proc_time)) + " (" + str(proc_time) + " seconds)")
    print("Log available at " + args.log.name)
    args.log.write("\nGenome slicing into oligos finished successfully at " + ctime() + "\n")
    args.log.write("Total time " + str(timedelta(seconds=proc_time)) + " (" + str(proc_time) + " seconds)\n")
