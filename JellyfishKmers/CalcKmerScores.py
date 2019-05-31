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

# Specify output filename
outputfile = "{filename}"

# Optional: Specify log filename (if none provided,
  defaults to output filename but with .log extension)
logfile = "{filename}"

# Calculate k-mer scores
exec(open("CalcKmerScores.py").read())
"""

import sys

# Custom exceptions
class NoInputError(NameError):
    pass
class NoOutputError(NameError):
    pass
class ScriptModeError(Exception):
    pass
class NoDictionaryError(Exception):
    pass


def Query(seq):
    level1 = seq[0:6]
    level2 = seq[0:12]
    count = counts[level1][level2][seq]
    return count

def CalcFromFasta(oligoname, outputname, dumpname, logname):
    # Setup IO
    source = open(oligoname, 'r')
    output = open(outputname, 'w')

    # Setup log file
    # If user specified log filename, use that
    if logname:
        pass
    # If printing output to screen, use generic log name
    elif outputname == "/dev/stdout":
        logname = "kmerscores.log"
    # Default: Same name as output file but with .log extension
    else:
        logname = outputname.rstrip(".fa") + ".log"
    log = open(logname, 'w')

    # Begin log file with context
    log.write("Log file for CalcKmerScores.py\n")
    log.write("Dictionary loaded from jellyfish dump file = " + dumpname + "\n")
    log.write("45-mer source file = " + oligoname + "\n")
    log.write("Output file of 45-mers and k-mer scores = " + outputname + "\n")

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

            # If 17-mer not found in dictionary, note in log and skip it
            except:
                log.write("No dictionary entry for " + seq + \
                " from source oligo " + header + "\n")
                continue

        output.write(header + " " + str(score) + "\n")
        output.write(oligo + "\n")
        header = source.readline().rstrip('\n')

    print("Finished writing k-mers and scores in fasta format to " + outputname)
    print("Log written to " + logname)

def CalcFromSam(samname, outputname, dumpname, logname):
    # Setup IO
    source = open(samname, 'r')
    output = open(outputname, 'w')

    # Setup log file
    # If user specified log filename, use that
    if logname:
        pass
    # If printing output to screen, use generic log name
    elif outputname == "/dev/stdout":
        logname = "kmerscores.log"
    # Default: Same name as output file but with .log extension
    else:
        logname = outputname.rstrip(".sam") + ".log"
    log = open(logname, 'w')

    # Begin log file with context
    log.write("Log file for CalcKmerScores.py\n")
    log.write("Dictionary loaded from jellyfish dump file = " + dumpname + "\n")
    log.write("Sam source file = " + samname + "\n")
    log.write("Output file of 45-mers and k-mer scores = " + outputname + "\n")

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

            # If 17-mer not found in dictionary, note in log and skip it
            except:
                log.write("No dictionary entry for " + seq + \
                " from source oligo " + splitline[0] + "\n")
                continue

        # Write line with k-mer score appended
        output.write(line.rstrip('\n') + "\tKS:i:" + str(score) + "\n")

        line = source.readline()

    print("Finished writing k-mers and scores in sam format to " + outputname)
    print("Log written to " + logname)


# ----------------main-------------------

try:
    # Make sure script is being run in interactive mode
    if sys.argv[0] == "CalcKmerScores.py":
        raise ScriptModeError

    # Make sure k-mer dictonary is loaded
    if not 'counts' in vars():
        raise NoDictionaryError

    if not 'outputfile' in vars():
        raise NoOutputError
    if not 'logfile' in vars():
        logfile = ""

    # Decide whether to use sam or fasta based on which variable is defined
    # Note: This program doesn't read anything from Jellyfish dump;
    # Passing Jellyfish dump filename to functions is solely for logging output
    if 'samfile' in vars():
        CalcFromSam(samfile, outputfile, dumpfile, logfile)
    elif 'oligofile' in vars():
        CalcFromFasta(oligofile, outputfile, dumpfile, logfile)
    else:
        raise NoInputError

except NoDictionaryError:
    print("Cannot calculate k-mer scores because k-mer dictionary not loaded.\n" + \
    "Please load the dictionary using the following commands:\n" + \
    \
    "dumpfile = \"{filename.fa}\"\n" + \
    "exec(open(\"LoadKmerDict.py\").read())\n" + \
    \
    "If you would like to read and write in fasta format, " + \
    "please set source filename as follows:\n" + \
    "oligofile = \"{filename.fa}\"\n" + \
    \
    "If you would like to read and write in sam format, " + \
    "please set source filename as follows:\n" + \
    "samfile = \"{filename.sam}\"\n" + \
    \
    "Please set output filename as follows:\n" + \
    "outputfile = \"{filename}\"\n")

except ScriptModeError:
    print("You ran CalcKmerScores.py in script mode.\n" + \
    "Please open the python shell and rerun in interactive mode " + \
    "using the following commands:\n" + \
    \
    "dumpfile = \"{filename.fa}\"\n" + \
    "exec(open(\"LoadKmerDict.py\").read())\n" + \
    \
    "If you would like to read and write in fasta format, " + \
    "please set source filename as follows:\n" + \
    "oligofile = \"{filename.fa}\"\n" + \
    \
    "If you would like to read and write in sam format, " + \
    "please set source filename as follows:\n" + \
    "samfile = \"{filename.sam}\"\n" + \
    \
    "Please set output filename as follows:\n" + \
    "outputfile = \"{filename}\"\n")

except NoOutputError:
    print("No output file specified.")
    print("Please set output filename as follows and try again:\n" + \
    "outputfile = \"{filename}\"")

except NoInputError:
    print("No input file specified.")
    print("If you would like to read and write in fasta format, " + \
    "please set filename as follows and try again:\n" + \
    "oligofile = \"{filename.fa}\"")
    print("If you would like to read and write in sam format, " + \
    "please set filename as follows and try again:\n" + \
    "samfile = \"{filename.sam}\"")

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
