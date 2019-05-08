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
try:
    import primer3
except ImportError:
    print("Primer3 not installed")
    exit(1)

# Filter out sequences with less than 70% homology
def bwa_filter(line, min_AS, max_XS):
    """
    :param line: line of sam file
    :param min_AS: min AS:i score, suggest probe length
    :param max_XS: max XS:i score, suggest probe length * homology
    :return: true to filter out line, false to keep
    """

    splitline = line.split('\t')
    try:
        AS_block = splitline[-2]
        XS_block = splitline[-1]
        assert AS_block[:5] == "AS:i:" and XS_block[:5] == "XS:i:" , "AS and XS not last two blocks"
    except AssertionError:
        for b in splitline:
            if b[:5] == "AS:i:":
                AS_block = b
            elif b[:5] == "XS:i:":
                XS_block = b

        try:
            assert AS_block[:5] == "AS:i:" and XS_block[:5] == "XS:i:" , "AS and XS not last two blocks"
            # print("Catastrophe avoided (last item was " + splitline[-1].rstrip('\n') + ")")
        except AssertionError:
            print("AS and XS not found in following line:")
            print(line)

    AS_score = int(AS_block[5:])
    XS_score = int(XS_block[5:])

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

    # Get sequence from line passed to function
    mapinfo = line.split('\t')
    sequence = mapinfo[9]

    # Calculate
    TM = primer3.calcTm(sequence)
    HTM = primer3.calcHairpinTm(sequence)

    # If melting temperature is too low, filter out
    if TM < min_TM:
        return True

    # If hairpin melting temperature is too high, filter out
    elif HTM > max_HTM:
        return True

    # If melting temperature and hairpin melting temperature
    # are too close together, filter out
    elif (TM - HTM) < min_diff_TM:
        return True

    else:
        return False

# Set up source and output files (returns tuple of file objects)
def setupio():
    # Read arguments from command line
    try:
        source_name = sys.argv[1]
        output_name = sys.argv[2]
        # rejected_name = sys.argv[3]

    # Use hard-codes file names if arguments not given
    except IndexError:
        source_name = "monkeywrench.sam"
        output_name = "dungeonmap_output2.sam"
        # rejected_name = "dungeonmap_rejects2.sam"
        # rejected_b_name = "dungeonmap_rejects_b2.sam"
        # rejected_p_name = "dungeonmap_rejects_p2.sam"

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
        if line[0] == '@':
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