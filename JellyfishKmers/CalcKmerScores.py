# 16 May 2019
# Lisa Malins
# CalcKmerScores.py

"""
Calculates k-mer scores for 45-mers using nested dictionary
of 17-mers and counts loaded in to memory by LoadKmerDict.py

From python shell:
dumpfile = "{filename.fa}"
exec(open("LoadKmerDict.py").read())
oligofile = "{filename.fa}"
exec(open("CalcKmerScores.py").read())
"""

def Query(seq):
    level1 = seq[0:6]
    level2 = seq[0:12]
    count = counts[level1][level2][seq]
    return count


#--------------- main -------------------
try:
    # Open source file of 45-mers
    oligofile = "fake45mers.fa"
    source = open(oligofile, 'r')

    line = source.readline().rstrip('\n')

    while line:
        # Error message for unreadable input
        assert line[0] == ">", \
        "\nUnable to read k-mers and scores due to unexpected input. " + \
        "Line was:\n" + line.rstrip('\n') + "\nfrom " + oligofile

        # Read 45-mer and calculate k-mer score
        line = source.readline().rstrip('\n')
        print("Calculating score for " + line)

        # Start with score of zero
        score = 0

        # Loop through 45-mer and query all 17-mers
        for i in range (0, 29):
            try:
                seq = line[i:i+17]
                count = Query(seq)
                score += int(count)

            # If 17-mer not found in dictionary, skip it
            except:
                print("No entry found for " + seq)
                continue

        print(line + " " + str(score))
        line = source.readline().rstrip('\n')


except FileNotFoundError:
    if oligofile:
        print("File " + oligofile + " not found")
    print("Please set filename as follows and try again:\n" + \
    "oligofile = \"{filename.fa}\"")
