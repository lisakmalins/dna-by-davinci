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

K-mer and header classes are in separate files.
Header class will need to be remade for different assemblies
    since header formats vary (see comments in class)
Header class lets this file know whether to get k-mers
    out of a section or skip it.

Usage:
python GetOligos.py {genome filename} {mer size} {step size} {output filename}

Example command:
python GetOligos.py agra_cadabra_genome.fa 45 3 agra_cadabra_45mers.fa
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

    # End program if end of file reached
    fo.close()
    sys.exit("Finished writing " + str(mer_size) + "-mers to " + output_name)





#-------------------main-----------------------

# Read arguments
try:
    source_name = sys.argv[1]
    mer_size = int(sys.argv[2])
    step_size = int(sys.argv[3])
    output_name = sys.argv[4]
    if mer_size <= 0 or step_size <= 0 or step_size > mer_size:
        raise ValueError
except IndexError:
    exit("Usage: python GetOligos.py {genome filename} {mer size} {step size} {output filename}")
except ValueError:
    exit("Usage: python GetOligos.py {genome filename} {mer size} {step size} {output filename}")


# Open source file in read-only mode
try:
    source = open(source_name, "r")
except FileNotFoundError:
    exit("File " + sys.argv[1] + " not found.")

# Set up file output
output = open(output_name, "w")

# Get length of file for progress output
print("Determining file length of " + source_name + "...")
source.seek(0,2)
filelength = float(source.tell())
source.seek(0)
percent = 10

print("Reading " + str(mer_size) + "-mers with step size of " + str(step_size) + \
" from " + source_name + " and writing to " + output_name)

# Create Kmer object
currentKmer = Kmer(source, mer_size, step_size)

# Loop through entire file
while not currentKmer.eof:
    # Output progress message
    if (source.tell() / filelength * 100) > percent:
        print("Read progress: " + str(percent) + "%")
        percent += 10

    # Find next header and start a new k-mer
    nextHeader = FindNextHeader(source)
    currentKmer.StartNewKMer(nextHeader.id)

    # Loop through sequence until next header
    while not currentKmer.eos:
        output.write(">" + currentKmer.id + "_" + str(currentKmer.index) + " \n" + str(currentKmer.seq) + "\n")
        currentKmer.GetNextKMer()

# Close file
source.close()
print("Finished writing " + str(mer_size) + "-mers to " + output_name)
