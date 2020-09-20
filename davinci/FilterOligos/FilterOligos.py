# 6 May 2019
# Lisa Malins
# FilterOligos.py

"""
Python program to filter hybridization oligos using bwa and primer3.
Based on bwa.py and primer3_filter.py from Chorus2 by zhangtaolab on GitHub.

Example command:
python FilterOligos.py -i unfiltered.sam -o filtered.sam

For full usage info, please see:
python FilterOligos.py --help
"""

import sys
from os import stat
try:
    import primer3
except ImportError:
    exit("primer3-py not installed")
from time import ctime
try:
    from time import process_time
except ImportError:
    from time import clock as process_time
from datetime import timedelta
import argparse

# Filter out sequences with less than 70% homology
def bwa_filter(line, min_AS, max_XS):
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
def primer3_filter(line, min_TM, max_HTM, min_diff_TM):
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
    parser = argparse.ArgumentParser(description="Filter oligos from SAM file based on BWA mapping statistics.\n")

    # I/O
    parser.add_argument("-i", "--in", dest="source", type=argparse.FileType('r'), help="input SAM filename", required=True)
    parser.add_argument("-o", "--output", type=argparse.FileType('w'), default="/dev/fd/1", help="args.output filename (default: standard out)")

    # BWA filtering
    parser.add_argument("--bwa-min-AS", dest="min_AS", type=int, default=45, help="minimum BWA alignment score (default: %(default)s)")
    parser.add_argument("--bwa-max-XS", dest="max_XS", type=int, default=31, help="maximum BWA suboptimal alignment score (default: %(default)s)")

    # Primer3 arguments
    parser.add_argument("--enable-primer3-filter", action="store_true", help="enable filtering by primer3 criteria")
    parser.add_argument("--min-TM", type=int, default=37, help="minimum melting temperature (default: %(default)s)")
    parser.add_argument("--max-HTM", type=int, default=35, help="maximum hairpin melting temperature (default: %(default)s)")
    parser.add_argument("--min-diff-TM", type=int, default=10, help="minimum difference between melting temperature and hairpin melting temperature (default: %(default)s)")

    # Other
    parser.add_argument("--write-rejects", action="store_true", help="write rejected oligos to separate output file")

    args = parser.parse_args()
    print(args)

    starttime = process_time()

    if (args.output.name == "/dev/fd/1"):
        log = open("/dev/null", 'w')
    else:
        log = open(args.output.name.rsplit('.', 1)[0] + ".log", 'w')
    log.write("Log file for FilterOligos.py")
    log.write("\nInput file to filter: " + args.source.name)
    log.write("\nFiltered args.output file: " + args.output.name)

    print("Filtering oligos from " + args.source.name + " and writing to " + args.output.name)

    if args.write_rejects:
        rejects = open(args.output.name.rsplit('.', 1)[0] + "_bwa_rejects.sam", 'w'), \
        open(args.output.name.rsplit('.', 1)[0] + "_primer3_rejects.sam", 'w')
        log.write("\nRejects written to: " + rejects[0].name + ", " + rejects[1].name)
        print("Writing rejects to " + rejects[0].name + ", " + rejects[1].name)

    # Echo arguments
    msg = "Filtering by parameters:"
    msg += "\nmin_AS = " + str(args.min_AS)
    msg += "\nmax_XS = " + str(args.max_XS)
    msg += "\nmin_TM = " + str(args.min_TM)
    msg += "\nmax_HTM = " + str(args.max_HTM)
    msg += "\nmin_diff_TM = " + str(args.min_diff_TM)

    print(msg)
    log.write("\n" + msg)

    # Setup status messages
    filelength = float(stat(args.source.name).st_size)
    percent = 10

    print("Filter beginning at " + ctime())
    log.write("\nFiltering began at " + ctime())

    # Loop through sam file
    while True:
        line = args.source.readline()
        if not line: break

        if (args.source.tell() / filelength * 100) >= percent:
            print("Progress: " + str(percent) + "% (" + ctime() + ")")
            percent += 10

        # Output all headers
        if line[0] == '@':
            args.output.write(line)
            continue

        # Discard lines that fail BWA filter
        elif bwa_filter(line, args.min_AS, args.max_XS):
            if args.write_rejects:
                rejects[0].write(line)
            continue

        # If primer3 filtering enabled, discard lines that fail primer3 filter
        elif args.enable_primer3_filter and primer3_filter(line, args.min_TM, args.max_HTM, args.min_diff_TM):
            if args.write_rejects:
                rejects[1].write(line)
            continue

        # Write lines that pass both filters
        else:
            args.output.write(line)


    # Close file
    args.source.close()

    endtime = process_time()
    proc_time = endtime - starttime

    msg = "Filtering completed successfully at " + ctime() + \
    "\nRun time: " + str(timedelta(seconds=proc_time)) + " (total seconds: " + str(proc_time) + ")"
    log.write("\n" + msg)
    print(msg)
    print("Filtered oligos written to " + args.output.name)
    print("Log written to " + log.name)
