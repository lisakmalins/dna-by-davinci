# 16 April 2019
# Lisa Malins
# GetOligos.py

"""
Script gets 45-mers out of .fa genome file.
Step size between 45-mers is 3 (values hardcoded)

Accepts as arguments the source filename, k-mer size, step size,
    and output filename.

Output is in fasta format:
    > [chromosome number] [index]
    [sequence]

K-mer and header classes are in separate files.
Header class will need to be remade for different assemblies
    since header formats vary (see comments in class)
Header class lets this file know whether to get k-mers
    out of a section or skip it.

Accepts as arguments the source filename, k-mer size, step size,
    and output filename.

Example command:
python GetOligos.py zmays_fake_genome.fa 45 3 test_output.fa
"""

import sys
from ZmaysB73Header import ZmaysHeader
from Kmer import Kmer


# Finds the next header and evaluates whether
#   we want to get k-mers out of that section.
# Ends program when no more relevant headers found.
def FindNextHeader(fo):
    # Read next line of file
    nextLine = fo.readline()

    # While there are more lines to read
    while nextLine:
        # If the next line is a header, create a header object and return to main
        if nextLine[0] == ">":
            header = ZmaysHeader(nextLine)

            # If the next line is a header that we care about,
            #   return the header and proceed to get k-mers from that sequence
            if header.EvalHeader():
                return header

            # If the next line is a header that we don't care about,
            #   skip it and keep looking
            else:
                # print("I don't care about ", header.id) #debug
                nextLine = fo.readline()
                continue

        # If the next line is not a header, read the next line and continue loop
        nextLine = fo.readline()

    # End program when end of file reached
    sys.exit(0)




#-------------------main-----------------------

# Read arguments
# (Comment out to use hard-coded argumemts)
source_name = sys.argv[1]
mer_size = int(sys.argv[2])
step_size = int(sys.argv[3])
output_name = sys.argv[4]

# Open source file in read-only mode
source = open(source_name, "r")
if source.closed:
    sys.exit("File open unsuccessful")

output = open(output_name, "w")

# Create Kmer object
currentKmer = Kmer(source, mer_size, step_size)

# Loop through entire file
while not currentKmer.eof:

    # Find next header and start a new k-mer
    nextHeader = FindNextHeader(source)
    currentKmer.StartNewKMer(nextHeader.id)

    # Loop through sequence until next header
    while not currentKmer.eos:
        output.write("> " + currentKmer.id + " " + str(currentKmer.index) + " \n" + str(currentKmer.seq) + "\n")
        currentKmer.GetNextKMer()

# Close file
source.close()
