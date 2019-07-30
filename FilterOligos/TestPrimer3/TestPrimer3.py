# 29 April 2019
# Lisa Malins
# TestPrimer3.py

"""
Independently runs primer3 on output from FilterOligos.py to evaluate
whether the oligos were kept or filtered for good reasons.

For each sequence, outputs melting temp, whether a hairpin was found,
the hairpin melting temp, and the difference between the two melting temps.

Includes smileys for the benefit of humans ;)

Example commands:
python3 TestPrimer3.py [sam file to evaluate] [where to output results]
python3 TestPrimer3.py test_output_p3only.sam /dev/stdout
python3 TestPrimer3.py test_rejects_p3only.sam /dev/stdout
"""

import sys
import primer3

# Convert true and false to :D and :(
def smiley(bool):
    if bool:
        return ":D"
    else:
        return ":("

# Set up source and output files (returns tuple of file objects)
def setupio():
    usage = "Usage: python3 TestPrimer3.py {sam file to evaluate}"
    try:
        source = open(sys.argv[1], 'r')
    except IndexError:
        exit(usage)
    except FileNotFoundError as e:
        exit("File " + e.filename + " not found.\n" + usage)

    try:
        output = open(sys.argv[2], 'w')
    except IndexError:
        output = open("/dev/stdout", 'w')

    return source, output


#-------------------main-----------------------

source, output = setupio()

# Iterate over lines in sam file passed to program
for line in source.readlines():
    # Skip headers
    if line[0] == '@':
        continue

    # Get sequence from sam file
    mapinfo = line.split('\t')
    sequence = mapinfo[9]

    # Calculate + evaluate melting temperature
    TM = primer3.calcTm(sequence)
    TM_smiley = smiley(TM > 37)

    # Calculate + evaluate hairpin
    hairpin = primer3.calcHairpin(sequence)
    hairpin_found = hairpin.structure_found
    HTM = hairpin.tm
    HTM_smiley = smiley(HTM < 35)

    # Calculate + evaluate melting temperature difference
    DTM = TM - HTM
    DTM_smiley = smiley(DTM > 10)

    # Evaluate the whole enchilada
    Overall_smiley = smiley((TM > 37) & (HTM < 35) & (DTM > 10))

    # Output results
    output.write(sequence)
    output.write("   TM= " + f"{TM:.01f} " + TM_smiley)
    output.write("   Hairpin= " + str(hairpin_found) + " " + f"{HTM:.01f} " + HTM_smiley)
    output.write("   DTM= " + f"{DTM:.01f} ".zfill(5) + DTM_smiley)
    output.write("   Overall = " + Overall_smiley)
    output.write('\n')

# Close file
source.close()
