# 16 May 2019
# Lisa Malins
# FakeFiles.py

"""
These 3 functions make synthetic test files for LoadKmerDict.py and CalcKmerScores.py.

fakedump() generates a fake file of 17-mers in the style of jellyfish dump.

fake45mers() generates a fake fasta-formatted file of 45-mers (alliteration
unintentional) in the style of Lisa's earlier masterpiece GetOligos.py.

fakesam() takes an existing sam file and swaps out 45-mers with fake ones.

Every 17-mer and 45-mer generated is a string of A's followed by a string of G's.
These files play nicely together because the 17-mers from fakedump() completely
cover the 45-mers from fake45mers() and fakesam().
"""

# Makes fake jellyfish dump file
def fakedump():
    output = open("fakedump.fa", 'w')

    # Start with 17-mer of A's with score of 1
    seq = "AAAAAAAAAAAAAAAAA"
    my_num = 1

    for x in range (0, 18):
        # Write score line
        output.write(">" + str(my_num) + "\n")
        # Write sequence
        output.write(seq + "\n")
        # Next sequence has 1 less A and 1 more G
        seq = seq[1:] + "G"
        # Next score is 1 greater
        my_num += 1

# Makes a fake file of 45-mers
def fake45mers():
    output = open("fake45mers.fa", 'w')
    # Start with 45-mer of A's with index 1 on chromosome 1
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    my_num = 1

    for x in range (0, 46):
        # Write header line
        output.write(">1_" + str(my_num) + "\n")
        # Write sequence line
        output.write(seq + "\n")
        # Next sequence has 1 less A and 1 more G
        seq = seq[1:] + "G"
        # Next index has step size of 3
        my_num += 3

# Takes existing sam file and swaps sequences for fake ones
def fakesam():
    originalsam = open("minimap.sam", 'r')
    output = open("fakemap.sam", 'w')

    # Start with 45-mer of A's
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

    for line in originalsam.readlines():
        # Print headers without touching them
        if line[0] == '@':
            output.write(line)
            continue

        # Replace real 45-mers with fake ones
        else:
            # Split line into components
            splitline = line.split('\t')

            # Reconstruct line with fake sequence
            fakeline = ""
            for i in range (0, 15):
                if i == 9:
                    fakeline += seq
                else:
                    fakeline += splitline[i]
                fakeline += '\t'
            output.write(fakeline.rstrip('\t'))

        # Next sequence has 2 less A's and 2 more G's
        seq = seq[2:] + "GG"
