# 6 May 2019
# Lisa Malins
# FilterOligos.py

"""
Python program to filter hybridization oligos using bwa and primer3.
Based on bwa.py and primer3_filter.py from Chorus2 by zhangtaolab on GitHub.

Example command:
python3 FilterOligos.py unfiltered.sam filtered.sam
"""

import sys
import re
try:
    import primer3
except ImportError:
    print("Primer3 not installed")
    exit(1)

# Find header lines and just print them
def isheader(line):
    # TODO can we not compile this every time?
    headerpat = re.compile('^@')
    if re.search(headerpat, line):
        return True
    else:
        return False

# Filter out sequences with less than 70% homology
def bwa_filter(line, min_AS, max_XS):
    """
    :param line: line of sam file
    :param min_AS: min AS:i score, suggest probe length
    :param max_XS: max XS:i score, suggest probe length * homology
    :return: true to filter out line, false to keep
    """

    # Generate patterns that will be used to search for
    # alignment scores (AS) and suboptimal alignment scores (XS)
    # TODO can we not compile this every time?
    AS_pat = re.compile('AS:i:(\d*)')
    XS_pat = re.compile('XS:i:(\d*)')

    # Search for an alignment score and a suboptimal alignment score
    # If found, save as match object (see re documentation)
    AS_match = re.search(AS_pat, line)
    XS_match = re.search(XS_pat, line)


    # If an alignment score and a suboptimal alignment score
    # are found, grab the scores and save for comparison to minimum and maximum.
    # If either score isn't found, return true (ignore this line)

    if AS_match:
        AS_score = int(AS_match.group(1))
    else:
        return True

    if XS_match:
        XS_score = int(XS_match.group(1))
    else:
        return True


    # If the alignment score is less than the minimum
    # or the suboptimal alignment score is greater than the maximum,
    # return true (filter out this sequence because we don't like it)
    if (AS_score < min_AS) or (XS_score >= max_XS):
        return True

    # Otherwise, if the scores are within range,
    # return false (do not filter)
    else:
        return False

# Filter out derpy oligos
def primer3_filter(line, min_TM=37, max_HTM=35, min_diff_TM=10):
    """
    :param sequence: sequence of oligo
    :param min_TM: minimum melting temp, default 37
    :param max_HTM: maximum hairpin melting temp, default 35
    :param min_diff_TM: minimum difference between melting temp and hairpin melting temp, default 10
    :return: true to filter out line, false to keep
    """
    primer3ft = False

    # Get sequence from line passed to function
    mapinfo = line.split('\t')
    sequence = mapinfo[9]

    # Calculate
    TM = primer3.calcTm(sequence)
    HTM = primer3.calcHairpinTm(sequence)

    # If melting temperature is too low, filter out
    if TM < min_TM:
        primer3ft = True

    # If hairpin melting temperature is too high, filter out
    if HTM > max_HTM:
        primer3ft = True

    # If melting temperature and hairpin melting temperature
    # are too close together, filter out
    if (TM - HTM) < min_diff_TM:
        primer3ft = True

    return primer3ft

# Set up source and output files (returns tuple of file objects)
def setupio():
    # Read arguments from command line
    try:
        source_name = sys.argv[1]
        output_name = sys.argv[2]
        # rejected_name = sys.argv[3]

    # Use hard-codes file names if arguments not given
    except IndexError:
        source_name = "dungeonmap.sam"
        output_name = "dungeonmap_output.sam"
        # rejected_name = "dungeonmap_rejects.sam"
        # rejected_b_name = "dungeonmap_rejects_b.sam"
        # rejected_p_name = "dungeonmap_rejects_p.sam"

    # Open source file in read-only mode
    source = open(source_name, 'r')

    # Set up file output
    output = open(output_name, 'w')
    # rejected = open(rejected_name, 'w')
    # rejected_b = open(rejected_b_name, 'w')
    # rejected_p = open(rejected_p_name, 'w')

    # return source, output, rejected, rejected_b, rejected_p
    return source, output

#-------------------main-----------------------

if __name__ == '__main__':

    # Set up source and output files
    # source, output, rejected, rejected_b, rejected_p = setupio()
    source, output = setupio()


# Iterate over lines in sam file passed to program
    for line in source.readlines():

        # Print headers without touching them
        if isheader(line):
            output.write(line)
            continue

        elif bwa_filter(line, 45, 31):
            # rejected.write(line.rstrip('\n') + " (bwa_filter)\n")
            # rejected_b.write(line)
            continue

        elif primer3_filter(line, 37, 35, 10):
            # rejected.write(line.rstrip('\n') + " (primer3_filter)\n")
            # rejected_p.write(line)
            continue

        else:
            output.write(line)


    # Close file
    source.close()