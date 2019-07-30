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
    exit("primer3-py not installed")

# Filter out sequences with less than 70% homology
def bwa_filter(line, min_AS=45, max_XS=31):
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

#-------------------main-----------------------

if __name__ == '__main__':

    # Set up source and output files
    usage = "Usage: python3 FilterOligos.py unfiltered.sam filtered.sam"

    # Read arguments from Snakefile
    try:
        source = open(snakemake.input[0], 'r')
        output = open(snakemake.output[0], 'w')
        write_rejected = snakemake.params.write_rejected
    # Or, read arguments from command line
    except NameError:
        try:
            source = open(sys.argv[1], 'r')
            output = open(sys.argv[2], 'w')
        except IndexError:
            exit(usage)
        except FileNotFoundError as e:
            exit("File " + e.filename + " not found.\n" + usage)

    # Default do not write rejects unless specified in Snakefile
    try:
        write_rejected = snakemake.params.write_rejected
    except NameError:
        write_rejected = False

    if write_rejected:
        rejects = open(output.name.rsplit('.', 1)[0] + "_bwa_rejects.sam", 'w'), \
        open(output.name.rsplit('.', 1)[0] + "_primer3_rejects.sam", 'w')
    else:
        rejects = False

    # Read parameters from Snakefile or use defaults
    try:
        min_AS = snakemake.params.bwa_min_AS
        max_XS = snakemake.params.bwa_max_XS
        min_TM = snakemake.params.primer3_min_TM
        max_HTM = snakemake.params.primer3_max_HTM
        min_diff_TM = snakemake.params.primer3_min_diff_TM
        print("Filtering by custom parameters:")
    except NameError:
        min_AS = 45
        max_XS = 31
        min_TM = 37
        max_HTM = 35
        min_diff_TM = 10
        print("Filtering by default parameters:")

    # Echo arguments to user
    print("min_AS = " + str(min_AS))
    print("max_XS = " + str(max_XS))
    print("min_TM = " + str(min_TM))
    print("max_HTM = " + str(max_HTM))
    print("min_diff_TM = " + str(min_diff_TM))

    print("Writing filtered sam to " + output.name)
    if write_rejected:
        print("Writing rejects to " + rejects[0].name + " " + rejects[1].name)

    # Iterate over lines in sam file passed to program
    for line in source.readlines():

        # Print headers without touching them
        if line[0] == '@':
            output.write(line)
            continue

        # Discard lines that do not pass BWA filter
        elif bwa_filter(line, min_AS, max_XS):
            if write_rejected:
                rejects[0].write(line)
            continue

        # Discard lines that do not pass primer3 filter
        elif primer3_filter(line, min_TM, max_HTM, min_diff_TM):
            if write_rejected:
                rejects[1].write(line)
            continue

        # Write lines that pass both filters
        else:
            output.write(line)


    # Close file
    source.close()
