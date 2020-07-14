import sys
import argparse
import re
try:
    import primer3
except ImportError:
    sys.stderr.write("primer3-py not installed")
    sys.exit(1)

"""
Returns reason for sequences that fail.
Returns False for good sequences.
"""
def primer3filter(seq, min_TM=37, max_HTM=35, min_diff_TM=10):

    # Calculate
    TM = primer3.calcTm(seq)
    HTM = primer3.calcHairpinTm(seq)

    # If melting temperature is too low, filter out
    if TM < min_TM:
        return "melting temp too low"

    # If hairpin melting temperature is too high, filter out
    elif HTM > max_HTM:
        return "hairpin melting temp too high"

    # If melting temperature and hairpin melting temperature
    # are too close together, filter out
    elif (TM - HTM) < min_diff_TM:
        return "difference between melting temp and hairpin melting temp too small"

    # If sequence will make a good probe, return false (do not filter)
    else:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter oligos from fasta. Check for homopolymers and primer3 criteria.\n")

    # I/O
    parser.add_argument("-i", "--in", dest="oligos", type=argparse.FileType('r'), help="input fasta filename", required=True)
    parser.add_argument("-o", "--output", type=argparse.FileType('w'), default="/dev/fd/1", help="output filename (default: standard out)")

    # Primer3 arguments
    parser.add_argument("--min-tm", type=int, default=37, help="minimum melting temperature (default: %(default)s)")
    parser.add_argument("--max-htm", type=int, default=35, help="maximum hairpin melting temperature (default: %(default)s)")
    parser.add_argument("--min-dtm", type=int, default=10, help="minimum difference between melting temperature and hairpin melting temperature (default: %(default)s)")

    parser.add_argument("--homopolymer-length", type=int, default=5, help="minimum length of homopolymer to filter out (default: %(default)s)")

    args = parser.parse_args()

    # Regex pattern which matches any sequence containing N
    # or a homopolymer of user-specified length or greater
    homopolymer = re.compile("N|A{{{n}}}|C{{{n}}}|G{{{n}}}|T{{{n}}}".format(n=args.homopolymer_length))

    # Loop through file
    linecount = 0
    while True:
        # Read and confirm header
        header = args.oligos.readline()
        if not header: break
        linecount += 1
        assert header[0] == ">", \
        "Oligo file {} not in recognized fasta format\nExpected fasta header on line {}, " \
        "instead found:\n{}\n".format(args.oligos.name, linecount, header)

        # Read sequence
        seq = args.oligos.readline()
        linecount += 1

        # Check for N's and homopolymers of 5 bases or more
        match = homopolymer.search(seq)
        if match:
            sys.stderr.write("Sequence {} failed homopolymer filter, sequence contains {}\n".format(header.strip(">\n"), match.group()))
            continue

        # Check for primer3 criteria
        p3filter = primer3filter(seq.rstrip(), args.min_tm, args.max_htm, args.min_dtm)
        if p3filter:
            sys.stderr.write("Sequence {} failed primer3 filter, reason: {}\n".format(header.lstrip(">").rstrip("\n"), p3filter))
            continue

        # Write sequence in fasta format if passes both filters
        args.output.write(header)
        args.output.write(seq)
