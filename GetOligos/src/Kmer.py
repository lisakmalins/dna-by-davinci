# 16 April 2019
# Lisa Malins
# Kmer.py

"""
Class that stores k-mers for GetOligos.py program.

Accepts .fa file to read k-mers from and uses counter
    to keep track of position.

Object remembers its k-mer sequence,
    what chromosome or other sequence it originates from,
    and its index within the chromosome.

Object also has flags for end of section
    and end of file to return to main.

Methods:
StartNewKMer() resets attributes and starts new 45-mer
    from next section in file.
GetNextKMer() forgets first 3 letters of sequence
    and fetches next 3 letters to get next k-mer (step size = 3)
    or returns to main if end of section or end of file reached.
"""

class Kmer:

    # Constructor
    def __init__(self, f):
        self.fo = f         # File to read from
        self.counter = 0    # Position in file (next letter to read)

        self.id = 0         # Chromosome ID
        self.index = 1      # Index within chromosome (refers to starting letter of seq)
        self.seq = ""       # K-mer sequence

        self.eos = False    # Has program reached the end of the section?
        self.eof = False    # Has program reached the end of the file?

    # To start new chromosome, reset attributes and get next 45-mer
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

        # Read next 45 letters, starting at stored counter position
        #   and ignoring non-alphabetical characters
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
            # If next character is >, set eos = true and return to main
            if nextChar:
                if nextChar.isalpha():
                    nextLetters += nextChar
                    self.index += 1

                elif nextChar == ">":
                    self.fo.seek(self.counter, 0)
                    self.eos = True
                    return

                self.counter += 1


            # If next character does not exist, exit function without appending to k-mer
            else:
                self.eof = True
                return

        # Forget first 3 letters of k-mer and append nextLetters
        self.seq = self.seq[3:] + nextLetters