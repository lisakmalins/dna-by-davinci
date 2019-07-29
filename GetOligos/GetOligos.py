# 16 April 2019
# Lisa Malins
# GetOligos.py

"""
Script gets k-mers out of fasta genome assembly.

Accepts as arguments the source filename, k-mer size, step size,
    and output filename.

Output is in fasta format:
    >[chromosome number]_[index]
    [sequence]

Usage:
python GetOligos.py {genome filename} {mer size} {step size} {output filename}

Example command:
python GetOligos.py agra_cadabra_genome.fa 45 3 agra_cadabra_45mers.fa
"""

import sys

class Kmer:

    # Constructor
    def __init__(self, f, k, s):
        self.fo = f         # File to read from
        self.counter = 0    # Position in file (next letter to read)

        self.k = k          # K-mer size
        self.s = s          # Step size (how far to advance on next k-mer)

        self.id = 0         # Chromosome ID
        self.index = 1      # Index within chromosome (refers to starting letter of seq)
        self.seq = ""       # K-mer sequence

        self.eos = False    # Has program reached the end of the section?
        self.eof = False    # Has program reached the end of the file?

    # To start new chromosome, reset attributes and get next k-mer
    def StartNewKMer(self, id):
        # Set chromosome id
        self.id = id

        # Set counter as current position of file object
        self.counter = self.fo.tell()

        # Reset index (position on chromosome) to 1
        self.index = 1

        # Reset end of section flag
        self.eos = False

        # Wipe stored k-mer sequence
        self.seq = ""

        # Read next k letters, starting at stored counter position
        #   and ignoring non-alphabetical characters
        while len(self.seq) < self.k:
            nextChar = self.fo.read(1)
            if nextChar.isalpha():
                self.seq += nextChar
            self.counter += 1

    # To get next k-mer, forget first s letters and append next s letters
    def GetNextKMer(self):
        nextLetters = ""
        self.fo.seek(self.counter, 0)

        # Fetch next s letters, ignoring spaces and newlines
        while len(nextLetters) < self.s:
            nextChar = self.fo.read(1)
            if nextChar:
                # If next character is a letter, add to nextLetters
                if nextChar.isalpha():
                    nextLetters += nextChar
                    self.index += 1
                # If next character is >, set eos = true and return to main
                elif nextChar == ">":
                    # Backspace file object's counter so header isn't skipped
                    self.fo.seek(self.counter, 0)
                    self.eos = True
                    return

                self.counter += 1


            # If next character does not exist, exit function without appending to k-mer
            else:
                self.eof = True
                return

        # Forget first s letters of k-mer and append nextLetters
        self.seq = self.seq[self.s:] + nextLetters


# Finds the next header and evaluates whether
#   we want to get k-mers out of that section.
# Ends program when no more relevant headers found.
def FindNextHeader(source, goodheaders):
    # Read next line of file
    line = source.readline()

    # While there are more lines to read
    while line:
        # If the next line is a header, create a header object and return to main
        if line[0] == ">":
            # If the next line is a header that we care about,
            #   return the header and proceed to get k-mers from that sequence
            if line in goodheaders:
                id = line[1:].split()[0]
                return id

            # If the next line is a header that we don't care about,
            #   skip it and keep looking
            else:
                # print("I don't care about ", header.id) #debug
                line = source.readline()
                continue

        # If the next line is not a header, read the next line and continue loop
        line = source.readline()

    # End program if end of file reached
    source.close()
    sys.exit("Finished writing " + str(mer_size) + "-mers to " + output.name)





#-------------------main-----------------------
usage = "Usage: python GetOligos.py {genome filename} {mer size} {step size} {output filename}"

# Read arguments
try:
    source = open(sys.argv[1], 'r')
except IndexError:
    exit(usage)
except FileNotFoundError as e:
    exit("File " + e.filename + " not found.\n" + usage)

try:
    mer_size = int(sys.argv[2])
    step_size = int(sys.argv[3])
    if mer_size <= 0 or step_size <= 0:
        raise ValueError()
    if step_size > mer_size:
        exit("Mer size must be greater than step size\n" + usage)
except IndexError:
    exit(usage)
except ValueError:
    exit("Please provide the mer-size and step-size as positive integers\n" + usage)

try:
    output = open(sys.argv[4], 'w')
except IndexError:
    output = open(source.name.split('.', 1)[0] + "_" + str(mer_size) + "mers.fa", 'w')

try:
    keep = open(source.name.rsplit('.',1)[0] + ".keep", 'r')
except FileNotFoundError as e:
    exit("File " + e.filename + " not found.\n" + \
    "See README on GitHub for instructions on how to specify which headers to keep for genome " + source.name)

# Read headers from keep file and verify they all have unique ID's
goodheaders = keep.readlines()
keep.close()

ids = set()
for h in goodheaders:
    id = h.lstrip('>').split()[0]
    if id in ids:
        message = "Error: sequence IDs must be unique. Fasta headers are improperly formatted. See main README for instructions on how to correct this.\n" #TODO

        message += "Headers known to program:\n"
        for h in goodheaders:
            message += h

        message += "IDs parsed so far:\n"
        for i in ids:
            message += i + "\n"

        message += "Duplicate: " + id + "\n"

        log.write(message)
        exit(message)

    else:
        ids.add(id)

# Get length of file for progress output
print("Determining file length of " + source.name + "...")
source.seek(0,2)
filelength = float(source.tell())
source.seek(0)
percent = 10

print("Reading " + str(mer_size) + "-mers with step size of " + str(step_size) + \
" from " + source.name + " and writing to " + output.name)

# Create Kmer object
currentKmer = Kmer(source, mer_size, step_size)

# Loop through entire file
while not currentKmer.eof:
    # Output progress message
    if (source.tell() / filelength * 100) > percent:
        print("Read progress: " + str(percent) + "%")
        percent += 10

    # Find next header and start a new k-mer
    id = FindNextHeader(source, goodheaders)
    currentKmer.StartNewKMer(id)

    # Loop through sequence until next header
    while not currentKmer.eos:
        output.write(">" + currentKmer.id + "_" + str(currentKmer.index) + " \n" + str(currentKmer.seq) + "\n")
        currentKmer.GetNextKMer()

# Close file
source.close()
print("Finished writing " + str(mer_size) + "-mers to " + output.name)
