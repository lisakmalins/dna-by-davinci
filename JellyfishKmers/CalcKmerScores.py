# 28 May 2019
# Lisa Malins
# CalcKmerScores.py

"""
Calculates k-mer scores for 45-mers using nested dictionary
of 17-mers and counts loaded in to memory by LoadKmerDict.py.
Reads and writes in either fasta format or sam format depending
on whether source filename is set with oligofile or samfile.

From python shell:
# Load k-mer scores from Jellyfish dump
dumpfile = "{filename.fa}"
exec(open("LoadKmerDict.py").read())

# Read k-mers from fasta or sam file (choose one)
oligofile = "{filename.fa}"
samfile = "{filename.sam}"

# Specify output filename (if none provided, defaults to standard out)
outputfile = "{filename}"

# Calculate k-mer scores
exec(open("CalcKmerScores.py").read())
"""


def Query(seq):
    level1 = seq[0:6]
    level2 = seq[0:12]
    count = counts[level1][level2][seq]
    return count

def CalcFromFasta(inputname, outputname):
    # Setup IO
    source = open(inputname, 'r')
    output = open(outputname, 'w')

    # Read 45-mers and calculate k-mer scores
    header = source.readline().rstrip('\n')
    while header:
        # Error message if line is not a fasta header
        assert header[0] == ">", \
        "\nUnable to read k-mers and scores due to unexpected input. " + \
        "Line was:\n" + line.rstrip('\n') + "\nfrom " + oligofile

        # Read 45-mers and calculate k-mer score
        oligo = source.readline().rstrip('\n')
        # print("Calculating score for " + line) #debug

        # Start with score of zero
        score = 0

        # Loop through 45-mer and query all 17-mers
        for i in range (0, 29):
            try:
                seq = oligo[i:i+17]
                count = Query(seq)
                score += int(count)

            # If 17-mer not found in dictionary, skip it
            except:
                print("No entry found for " + seq)
                continue

        output.write(header + " " + str(score) + "\n")
        output.write(oligo + "\n")
        header = source.readline().rstrip('\n')

    print("Finished writing k-mers and scores in fasta format to " + outputname)

def CalcFromSam(inputname, outputname):
    # Setup IO
    source = open(inputname, 'r')
    output = open(outputname, 'w')

    # Read 45-mers and calculate k-mer scores
    line = source.readline()
    while line:
        # Print headers without touching them
        if line[0] == '@':
            output.write(line)
            line = source.readline()
            continue

        # Split line into components and grab oligo
        splitline = line.split('\t')
        oligo = splitline[9]

        # Start with score of zero
        score = 0

        # Loop through 45-mer and query all 17-mers
        for i in range (0, 29):
            try:
                seq = oligo[i:i+17]
                count = Query(seq)
                score += int(count)

            # If 17-mer not found in dictionary, skip it
            except:
                print("No entry found for " + seq)
                continue

        # Write line with k-mer score appended
        output.write(line.rstrip('\n') + "\tKS:i:" + str(score) + "\n")

        line = source.readline()

    print("Finished writing k-mers and scores in sam format to " + outputname)


# ----------------main-------------------

try:
    # If output file not specified, print to screen
    if not 'outputfile' in vars():
        outputfile = "/dev/stdout"

    # Decide whether to use sam or fasta based on which variable is defined
    if 'samfile' in vars():
        CalcFromSam(samfile, outputfile)
    elif 'oligofile' in vars():
        CalcFromFasta(oligofile, outputfile)
    else:
        raise NameError


except NameError:
    print("No input and/or output file specified.")
    print("If you would like to read and write in fasta format, " + \
    "please set filename as follows and try again:\n" + \
    "oligofile = \"{filename.fa}\"")
    print("If you would like to read and write in sam format, " + \
    "please set filename as follows and try again:\n" + \
    "samfile = \"{filename.sam}\"")
    # print("Please set output filename as follows and try again:\n" + \
    # "outputfile = \"{filename}\"")
except FileNotFoundError:
    if 'oligofile' in vars():
        print("File " + oligofile + " not found\n" + \
        "Please set filename as follows and try again:\n" + \
        "oligofile = \"{filename.fa}\"")
    elif 'samfile' in vars():
        print("File " + samfile + " not found\n" + \
        "Please set filename as follows and try again:\n" + \
        "samfile = \"{filename.sam}\"")
    # print("Please set output filename as follows and try again:\n" + \
    # "outputfile = \"{filename}\"")
