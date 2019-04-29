# 28 April 2019
# Lisa Malins
# CheckForRepeats.py

"""
Checks for repeats in a sam file.
Accepts source file as command-line argument.
Prints names and locations of repeated sequences or "No repeats found".
"""

import sys
sys.path.append('../')
from FilterOligos import isheader

# Read arguments (comment out to hard-code values)
source_name = sys.argv[1]

# # Hard-code file name (comment out to accept filename as argument)
# source_name = "CheckForRepeats/positive_control.sam"
# source_name = "tinymap.sam"

# Open source file in read-only mode
source = open(source_name, 'r')
if source.closed:
    sys.exit("File open unsuccessful")

# Initialize variables used in loop
last_id = ""
last_seq = ""
next_id = ""
next_seq = ""
repeats_found = 0

for line in source.readlines():
    # Skip headers
    if isheader(line):
        continue

    # Grab sequence ID
    mapinfo = line.split('\t')
    id = mapinfo[0]
    seq = mapinfo[2]

    # Transfer last id and sequence read
    last_id = next_id
    last_seq = next_seq

    # Read next id and sequence
    next_id = id
    next_seq = seq

    if last_id == next_id:
        print("ID " + next_id + " is repeated in sequences " + last_seq + " and " + next_seq)
        repeats_found += 1


if not repeats_found:
    print("No repeats found in " + source_name)
else:
    print(str(repeats_found) + " repeats found in " + source_name)