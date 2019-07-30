# 6 May 2019
# Lisa Malins
# FilterOligos.py

"""
Python program to filter hybridization oligos using bwa and primer3.
Based on bwa.py and primer3_filter.py from Chorus2 by zhangtaolab on GitHub.

Filtering ranges can be configured with Snakemake.

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
    :param min_AS: minimum alignment score (sam tag AS:i:), ideally probe length
    :param max_XS: maximum suboptimal alignment score (sam tag XS:i), ideally probe length * homology
    :return: true to filter out line, false to keep
    """

    fields = line.split('\t')

    # AS and XS usually last two fields but not always
    try:
        AS_field = fields[-2]
        XS_field = fields[-1]
        assert AS_field[:5] == "AS:i:" and XS_field[:5] == "XS:i:" , "AS and XS not last two blocks"

    except AssertionError:
        # Check all fields if necessary
        for f in fields:
            if f[:5] == "AS:i:":
                AS_field = f
            elif f[:5] == "XS:i:":
                XS_field = f
        try:
            assert AS_field[:5] == "AS:i:" and XS_field[:5] == "XS:i:"
        except AssertionError:
            print("AS and XS not found in following line:")
            print(line)
            return True

    try:
        AS_score = int(AS_field[5:])
        XS_score = int(XS_field[5:])
    except ValueError:
        print("Could not parse integers from AS and XS in following line:")
        print(line, AS_field, XS_field)
        return True

    # If the alignment score < minimum or suboptimal alignment score > maximum,
    # return true (filter this sequence because it maps poorly)
    if (AS_score < min_AS) or (XS_score >= max_XS):
        return True

    # If the scores are within range, return false (do not filter)
    else:
        return False

# Filter out oligos that would behave unexpectedly as probes
def primer3_filter(line, min_TM=37, max_HTM=35, min_diff_TM=10):
    """
    :param sequence: sequence of oligo
    :param min_TM: minimum melting temp, default 37
    :param max_HTM: maximum hairpin melting temp, default 35
    :param min_diff_TM: minimum difference between melting temp and hairpin melting temp, default 10
    :return: true to filter out line, false to keep
    """

    # Get sequence from line passed to function
    sequence = line.split('\t')[9]

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

    # If sequence will make a good probe, return false (do not filter)
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


    # Loop through sam file
    for line in source.readlines():

        # Output all headers
        if line[0] == '@':
            output.write(line)
            continue

        # Discard lines that fail BWA filter
        elif bwa_filter(line, min_AS, max_XS):
            if write_rejected:
                rejects[0].write(line)
            continue

        # Discard lines that fail primer3 filter
        elif primer3_filter(line, min_TM, max_HTM, min_diff_TM):
            if write_rejected:
                rejects[1].write(line)
            continue

        # Write lines that pass both filters
        else:
            output.write(line)


    # Close file
    source.close()
