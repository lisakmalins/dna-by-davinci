# 23 April 2019
# Lisa Malins
# FilterOligos.py

"""
Python program to filter hybridization oligos using bwa and primer3.
Based on bwa.py and primer3_filter.py from Chorus2 by zhangtaolab on GitHub.
"""

import sys
import re
try:
    import primer3
except ImportError:
    pass

# Find header lines and just print them
def isheader(line):
    # TODO can we not compile this every time?
    headerpat = re.compile('^@')
    if re.search(headerpat, line):
        return True
    else:
        return False

# Filter out sequences with less than 70% homology
def bwafilter(line, min_AS, max_XS):
    """
    :param line: line of sam file
    :param min_AS: min AS:i score, suggest probe length
    :param max_XS: max XS:i score, suggest probe length * homology
    :return: true to filter out line, false to keep
    """

    # Generate patterns that will be used to search for
    # alignment scores (AS) and suboptimal alignment scores (XS)
    # TODO can we not compile this every time?
    AS_pat = re.compile('AS:i:(\d.)')
    XS_pat = re.compile('XS:i:(\d.)')

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

    # print("AS = " + str(AS_score) + " min_AS = " + str(min_AS) + " XS = " + str(XS_score) + " max_XS =" + str(max_XS))


    # If the alignment score is less than the minimum
    # or the suboptimal alignment score is greater than the maximum,
    # return true (filter out this sequence because we don't like it)
    if (AS_score <= min_AS) & (XS_score > max_XS):

        # mapinfo = line.split('\t')
        #
        # seqlist.append(mapinfo[9])

        # print("AS score of " + str(AS_score) + " is greater than or equal to " + str(min_AS))
        # print("XS score of " + str(XS_score) + " is less than" + str(max_XS))

        return True

    # Otherwise, if the scores are within range,
    # return false (do not filter)
    else:
        return False


# Screens derpy oligos by making sure that: melting temperature is above the minimum,
# hairpin melting temperature is below the maximum, and the difference between the two is sufficient.
# Returns boolean: true = filter out, false = keep.
def primer3_filter(sequence, mintm=37, maxhtm=35, dtm=10):

    primer3ft = False

    tm = primer3.calcTm(sequence)

    htm = primer3.calcHairpinTm(sequence)

    ## If melting temperature is too low, filter out
    if tm < mintm:
        primer3ft = True

    ## If hairpin melting temperature is too high, filter out
    if htm > maxhtm:
        primer3ft = True

    ## If melting temperature and hairpin melting temperature
    ## are too close together, filter out
    if (tm-htm) < dtm:
        primer3ft = True

    # print(sequence, tm, htm, dtm)

    return primer3ft


#-------------------main-----------------------

if __name__ == '__main__':

    # # Read arguments (comment out to hard-code values)
    # source_name = sys.argv[1]
    # output_name = sys.argv[2]

    # Hard-code file names (comment out to accept as argument)
    source_name = "tinymap.sam"
    output_name = "test_output.sam"

    # Open source file in read-only mode
    source = open(source_name, 'r')
    if source.closed:
        sys.exit("File open unsuccessful")

    # Set up file output
    output = open(output_name, 'w')

    ## Iterate over lines in sam file passed to program
    for line in source.readlines():

        # Print headers without touching them
        if isheader(line):
            print(line.rstrip('\n'))
            continue

        elif bwafilter(line, 45, 31):
            print("We don't like this sequence:")
            print(line.rstrip('\n'))

        else:
            print("We like this sequence:")
            print(line.rstrip('\n'))


    # Close file
    source.close()