# 3 June 2019
# Lisa Mains
# GetSizeOfDict.py
"""
Recursive function that approximates memory footprint of nested dictionary.
Adds sizes of dictionaries, their keys, and their values.

Usage:
from GetSizeOfDict import sizeofdict
sizeofdict( {dictionary name} )
"""

import sys

def sizeofdict(d, sum=0, verbose=False):
    # Execute this block for innermost values only
    # Add size of final value and stop recursion
    if not isinstance(d, dict):
        # Add size of final value
        sum += sys.getsizeof(d)
        if verbose:
            print("Size of value", d, "is", + sys.getsizeof(d))
        return sum
    # Execute this block for dictionaries
    # For each item, add size of key and recur on members
    for item in d:
        # Add size of string key
        sum += sys.getsizeof(item)
        if verbose:
            print("Size of key", item, "is", str(sys.getsizeof(item)))
        # Recur on members of dict
        sum = sizeofdict(d[item], sum, verbose)

    # Add size of dictionary
    sum += sys.getsizeof(d)
    if verbose:
        print("Size of ", type(d), " is ", str(sys.getsizeof(d)))
    # print("Current total is ", sum)
    return sum

# Example usage
if __name__ == '__main__':
    TrekCharacters = {'TOS': {'Kirk': 'Captain', 'Spock': 'Commander'}, \
     'DS9': {'Sisko': 'Captain', 'Dax': 'Lieutenant'}}
    print(sizeofdict(TrekCharacters, verbose=True))
