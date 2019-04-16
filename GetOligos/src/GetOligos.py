# 16 April 2019
# Lisa Malins
# GetOligos.py

"""
Script gets 45-mers out of .fa genome file.
Step size between 45-mers is 3 (values hardcoded)

Program removes newlines from 45-mers.
Output is [chromosome number],[index],[sequence].

Includes easy comment-out swap to switch between
    accepting filename as command-line argument
    and hard-coding filename to run quickly from IDE
    or from command-line without arguments.

K-mer and header classes are in separate files.
Header class will need to be remade for different assemblies
    since header formats vary (see comments in class)
Header class lets this file know whether to get k-mers
    out of a section or skip it.
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

# Read file name as argument
# (Comment out to use hard-coded filename)
filename = sys.argv[1]

# Hard-code filename
# (Comment out to accept file as command-line argument)
# filename = "zmays_fake_genome.fa"

# Open file in read-only mode
fo = open(filename, "r")
if fo.closed:
    sys.exit("File open unsuccessful")

# Create Kmer object
currentKmer = Kmer(fo)

# Loop through entire file
while not currentKmer.eof:

    # Find next header and start a new k-mer
    nextHeader = FindNextHeader(fo)
    currentKmer.StartNewKMer(nextHeader.id)

    # Loop through sequence until next header
    while not currentKmer.eos:
        print(currentKmer.id, ",", currentKmer.index, ",", currentKmer.seq)
        currentKmer.GetNextKMer()

# Close file
fo.close()