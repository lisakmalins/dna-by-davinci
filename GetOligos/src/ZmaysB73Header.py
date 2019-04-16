# 16 April 2019
# Lisa Malins
# ZMaysB73Header.py

"""
Class that stores information extracted from headers of
Zea mays ssp mays cv B73 Reference Genome
from <https://www.maizegdb.org/assembly>.

Since header format varies between genome assemblies,
a different version of this file will need to be written
for another assembly with a different header format.

Necessary components for rewrite:
    id -- string that holds identification for a given sequence
        (i.e., chromosome number, whether it's from a plastid or
        mitochondrion, or contig identification)
    EvalHeader() -- returns True or False whether main should
        get k-mers out of this sequence or not.
"""


class ZmaysHeader:

    # Constructor -- accepts header text and extracts info
    def __init__(self, line):
        self.line = line                      # entire line

        self.words = line.split()
        self.type = self.words[1]             # dna:chromosome or dna:contig
        self.id = self.words[0][1:]           # e.g., chromosome number or Pt, Mt, B73V4_ctg150

        self.words = self.words[2].split(":")
        self.length = self.words[4]           # length of sequence

    # Evaluate whether we want to read k-mers from a section based on its header
    # If the header indicates section is a normal, numbered chromosome, returns true
    # Function needs to be rewritten if different assembly has different header format
    def EvalHeader(self):
        if self.type == "dna:chromosome" and self.id.isdigit():
            return True
        else:
            return False

    # Debugging function: Prints all attributes
    def PrintAll(self):
        print("Type: ", self.type)
        print("ID: ", self.id)
        print("Length: ", self.length)
        print("Do we care: ", self.EvalHeader(), "\n")
