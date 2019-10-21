# 29 July 2019
# Lisa Malins
# GetOligos.py

"""
Script slices fasta genome assembly into overlapping oligos.

Accepts as arguments the source filename, k-mer size, step size,
    output filename (optional) and log filename (optional).

Usage:
python GetOligos.py {genome filename} {keep sequences filename} {mer size} {step size} {optional: output filename} {optional: log filename}

Example command:
python GetOligos.py agra_cadabra_genome.fa 45 3 agra_cadabra_45mers.fa
"""

import sys
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


# Finds headers and returns "ID" next header that that matters
def NextHeader(source, log, good_seqs, line=" "):
    while True:
        # If the next line is a header
        if line[0] == ">":
            id = line[1:].split()[0]
            if id in good_seqs:
                return id
            else:
                log.write("Ignored:\n" + line)
                # print("I don't care about ", line.rstrip()) #debug

        line = source.readline()

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
            print(next_letter + " <-- wtf")
            exit("I literally can't even right now") #debug
    return addition

def SetupIO():
    usage = "Usage: python GetOligos.py {genome filename} {keep sequences filename} {mer size} {step size} {optional: output name} {optional: log name}"

    # Read arguments
    try:
        source = open(sys.argv[1], 'r')
    except IndexError:
        exit(usage)
    except FileNotFoundError as e:
        exit("Genome file " + e.filename + " not found.\n" + usage)

    try:
        keep = open(sys.argv[2], 'r')
    except FileNotFoundError as e:
        keep_usage = "This should be a text file of the sequences you want to get oligos from, one per line.\n"

        # Attempt to catch user forgetting keep file
        if sys.argv[2].isdigit():
            sys.stderr.write("I'm guessing you forgot to include a file of sequences to keep.\n{}{}\
            ".format(keep_usage, usage))
            exit(1)

        # User specified a keep file but doesn't exist
        sys.stderr.write("Sequences to keep file " + e.filename + " not found.\n{}{}\
        ".format(keep_usage, usage))
        exit(1)

    try:
        mer_size = int(sys.argv[3])
        step_size = int(sys.argv[4])
        if mer_size <= 0 or step_size <= 0:
            raise ValueError()
        if step_size > mer_size:
            exit("Mer size must be greater than step size\n" + usage)
    except IndexError:
        exit(usage)
    except ValueError:
        exit("Please provide the mer-size and step-size as positive integers\n" + usage)

    try:
        output = open(sys.argv[5], 'w')
    except IndexError:
        output = open(source.name.rsplit('.', 1)[0] + "_" + str(mer_size) + "mers.fa", 'w')

    try:
        log = open(sys.argv[6], 'w')
    except IndexError:
        log = open(output.name.rsplit('.', 1)[0] + ".log", 'w')

    return source, mer_size, step_size, output, log, keep


#-------------------main-----------------------

source, mer_size, step_size, output, log, keep = SetupIO()

# Read headers from keep file and verify they all have unique ID's
good_seqs = list(map(lambda x: str(x).strip(), keep.readlines()))
keep.close()

log.write("Log file for GetOligos.py\n")
log.write("Genome file: " + source.name + "\n")
log.write("Oligo size: " + str(mer_size))
log.write(" Step size: " + str(step_size) + "\n")
log.write("Sequences to keep read from: " + keep.name + "\n")
log.write("\n".join(good_seqs))
log.write("\n")
log.write("Oligos written to: " + output.name + "\n\n")


# Get length of file for progress output
filelength = float(stat(source.name).st_size)
percent = 10

print("Reading " + str(mer_size) + "-mers with step size of " + str(step_size) + \
" from " + source.name + " and writing to " + output.name)
print("Logging to: " + log.name)
print("Genome slicing into oligos beginning at " + ctime())
log.write("\nGenome slicing into oligos beginning at " + ctime() + "\n\n")
time0 = process_time()

# Create Kmer object
kmer = Kmer()
line = " "

try:
    while True:
        # Progress messages
        if (source.tell() / filelength * 100 >= percent):
            print("Read progress : " + str(percent) + "%")
            percent += 10

        # Get header, start new k-mer
        id = NextHeader(source, log, good_seqs, line)
        kmer.StartNew(id, ReadChars(source, 45))

        # Print k-mers until header encountered
        while True:
            output.write(">" + str(id) + "_" + str(kmer.index) + "\n")
            output.write(kmer.seq + "\n")
            try:
                kmer.Advance(ReadChars(source, 3))
            # Catch header line and carry to next main loop iteration
            except HeaderException as e:
                line = e.args[0]
                break

# Aaaaaaand stick the landing
except (IndexError, EOFError) as e:
    source.close()

    proc_time = process_time() - time0

    print("Finished writing " + str(mer_size) + "-mers to " + output.name)
    print("Program finished successfully at " + ctime())
    print("Total time " + str(timedelta(seconds=proc_time)) + " (" + str(proc_time) + " seconds)")
    print("Progress messages may not have reached 100% because of ignored sequences.")
    print("Log available at " + log.name)
    log.write("\nGenome slicing into oligos finished successfully at " + ctime() + "\n")
    log.write("Total time " + str(timedelta(seconds=proc_time)) + " (" + str(proc_time) + " seconds)\n")
