# 2 April 2019
# Lisa Malins
# GetOligos.py

"""
Simple script that gets 45-mers out of .fa file without chromosome headers.
Step size between 45-mers is 3.
Program removes newlines from 45-mers.
Output is 45-mer sequence only.
Accepts filename as command-line argument.
Algorithm is slow because the program forgets the 45-mer every time
    and fetches a new one using the counter,
    re-checking entire 45-mer for newlines.
Newer version uses faster algorithm.
"""

import sys

# Print all 45-mers
def Print45Mers(f):
    # Open file
    fo = open(f, "r") #open file in read-only mode
    if fo.closed:
        print("File open unsuccessful")

    # Get size of file
    fo.seek(0,2) # move to end of file
    size = fo.tell()
    #print("Size of file is : ", size)

    # Print all 45-mers in file
    counter = 0 # keeps track of file index
    while counter < size - 45:
        fo.seek(counter,0) # set seek position to counter
        str = fo.read(45) # set str to the next 45-mer

        #correct strings with linebreaks
        if str.find("\n") != -1:
            str = str.replace("\n", "") # remove linebreak
            str = str + fo.read(1) # append next char

        print(str)
        counter = counter + 3 # increment counter by 3

    # Close file
    fo.close()

filename = sys.argv[1]
Print45Mers(filename)
