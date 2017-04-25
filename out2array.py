# -----------------------------------------------------------------
# out2array.py
#
# This module allows us to format the files output from mixed_ising
# into a style that is easy to feed in to mathematica directly.
# -----------------------------------------------------------------

import sys
import re
import os

def read_lines(lines, points):
    inclusion = re.compile("is not excluded")   # regex compile for mixed_ising output
    number_data = re.compile("[\d]+.[\d]+")     # general format for a double outputted
    for line in lines:
        if inclusion.search(line) is not None:  # check that this point is included
            data = number_data.findall(line)    # find the 2 or 3 coordinates of this point
            points.append(map(lambda x: float(x), data)) # add it as an array of floats

if __name__ == "__main__":

    points = []
    # If one argument is supplied, it'll search through the out files
    # for all those matching that argument
    if len(sys.argv) > 1:
        for filename in os.listdir("out_files"):
            if (sys.argv[1]+'_') in filename:
                f_out = open("out_files/{}".format(filename), 'r')
                read_lines(f_out, points)
    # Otherwise we read from the file that is fed into the program
    else:
        read_lines(sys.stdin.readlines(), points)

    output = "{}".format(points)        # Turns array of points to a string
    output = output.replace("[", "{")   # Replaces the square brackets
    output = output.replace("]", "}")   # to give Mathematica-friendly format

    print output
