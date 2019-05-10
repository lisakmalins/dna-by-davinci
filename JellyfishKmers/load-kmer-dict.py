# From python shell:
# import sys
# sys.argv = ['filename.fa']
# exec(open("load-kmer-dict.py").read())
# sys.getsizeof(scores)

try:
    sourcename = sys.argv[0]
    source = open(sourcename, 'r')
    print("Reading kmer scores from file " + sys.argv[0] + "...")
    line = source.readline()
    scores = dict()
    while line:
        assert (line[0] == ">")
        score = line[1:].rstrip('\n')
        line = source.readline()
        seq = line.rstrip('\n')
        scores[seq] = score
        line = source.readline()
    print(str(len(scores)) + " kmers and scores read from file " + sys.argv[0])
    print("\tEnter scores to print all entries in dictionary")
    print("\tEnter sys.getsizeof(scores) to find size of dictionary")
    print("\tEnter scores['sequence'] to see score for a particular sequence")

except NameError:
    print("Please import sys, set filename as sys.argv[0], and try again")
    print("Make sure you are running this within the python shell or there is no point")
except FileNotFoundError:
    if sys.argv[0]:
        print("File " + sys.argv[0] + " not found")
    print("Please set filename as sys.argv[0] and try again")
except AssertionError:
    print("Unexpected input. Line was: " + line.rstrip('\n'))
