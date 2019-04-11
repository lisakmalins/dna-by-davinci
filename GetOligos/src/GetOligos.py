# 8 April 2019
# Lisa Malins
# GetOligos.py

"""
Script gets 45-mers out of .fa chromosome chunk (without headers).
Step size between 45-mers is 3.
Program removes newlines from 45-mers.
Output is [chromosome number] [index] [sequence].
Includes easy comment-out swap to switch between
    accepting filename as command-line argument
    and hard-coding filename to run quickly from IDE (IntelliJ).
New algorithm is faster than previous versions because instead of
    forgetting each 45-mer and starting a new one every time using
    a counter, the program discards the first 3 letters only and
    appends the next 3 letters in the sequence.
"""

import sys

class Kmer:

    # Constructor
    def __init__(self, f):
        self.fo = f         # File to read from
        self.index = 1      # Index within chromosome (refers to starting letter of seq)
        self.counter = 0    # Position in file (next letter to read)
        self.seq = ""       # K-mer sequence
        self.chrom = 0      # Chromosome number
        self.eof = False    # Has program reached the end of the file?

    # To start new chromosome
    def StartNewKMer(self):
        # Increment chromosome number
        self.chrom += 1

        # Reset index (position on chromosome) to 1
        self.index = 1

        # Read next 45 letters, starting at stored counter position
        #   and ignoring non-alphabetical characters
        self.fo.seek(self.counter, 0)
        while len(self.seq) < 45:
            nextChar = self.fo.read(1)
            if nextChar.isalpha():
                self.seq += nextChar
            self.counter += 1


    # To get next k-mer, forget first 3 letters and append next 3 letters
    def GetNextKMer(self):
        nextLetters = ""
        self.fo.seek(self.counter, 0)

        # Fetch next 3 letters, ignoring spaces and newlines
        while len(nextLetters) < 3:
            nextChar = self.fo.read(1)

            # If next character exists and is a letter, add to nextLetters
            # Python reads an empty string as false, anything else is true
            if nextChar:
                if nextChar.isalpha():
                    nextLetters += nextChar
                    self.index += 1
                self.counter += 1

            # If next character does not exist, exit function without appending to k-mer
            else:
                self.eof = True
                return

        # Forget first 3 letters and append nextLetters
        self.seq = self.seq[3:] + nextLetters




#-------------------main-----------------------


# Read file name as argument
# (Comment out to use hard-coded filename)
filename = sys.argv[1]

# Hard-code filename to run quickly in IntelliJ
# (Comment out to accept file as command-line argument)
# filename = "Zea_mays_testsample2.fa"

# Open file
fo = open(filename, "r")  # open file in read-only mode
if fo.closed:
    sys.exit("File open unsuccessful")

# # Get size of file
# fo.seek(0,2) # move to end of file
# size = fo.tell()
# #print("Size of file is : ", size)

# Start new chromosome
currentKmer = Kmer(fo)
currentKmer.StartNewKMer()

# Print all k-mers
# Print executes before GetNextKMer() in order to print
#   1st k-mer from StartNewKmer() but not print final
#   incomplete k-mer if number of letters is not a multiple of 3
while not currentKmer.eof:
    # Print [chromosome number] [index] [sequence] separated by tabs
    print(currentKmer.chrom, "	", currentKmer.index, "	", currentKmer.seq)
    currentKmer.GetNextKMer()

# Close file
fo.close()