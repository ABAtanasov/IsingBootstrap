# -----------------------------------------------------------------
# out2array.py
#
# This module allows us to format the files output from mixed_ising
# into a style that is easy to feed in to mathematica directly.
# -----------------------------------------------------------------

import sys
import re
import os
from point_generator import generate_from_file

def read_lines(lines, points, form="is not excluded"):
    inclusion = re.compile("(\([-?\d.\d, ]+\)) " + form)
    number_data = re.compile("-?[\d]+.[\d]+")     # general format for a double outputted
    for line in lines:
        inclusionstring = inclusion.search(line)
        if inclusionstring is not None:                              # check that this point is included
            data = number_data.findall(inclusionstring.groups()[0])  # get the 2 or 3 coordinates of the point
            points.append(map(lambda x: float(x), data)) # add it as an array of floats

if __name__ == "__main__":

    points = []
    if len(sys.argv) > 1 and sys.argv[1] == "envelope":
        envelope = True
    else:
        envelope = False
    # If one argument is supplied, it'll search through the out files
    # for all those matching that argument
    if len(sys.argv) > 1:
        form = "is not excluded"
        if envelope:
            form = "is excluded"

        for filename in os.listdir("out_files"):
            if (sys.argv[-1]+'_') in filename:
                f_out = open("out_files/{}".format(filename), 'r')
                read_lines(f_out, points, form=form)
                f_out.close()

    # Otherwise we read from the file that is fed into the program
    else:
        read_lines(sys.stdin.readlines(), points)


    if envelope:
        with open("in_files/{}.pts".format(sys.argv[-2]), 'r') as f_in:
            all_points = generate_from_file(f_in=f_in)
            points = [point for point in all_points if point not in points]


    output = "{}".format(points)        # Turns array of points to a string
    output = output.replace("[", "{")   # Replaces the square brackets
    output = output.replace("]", "}")   # to give Mathematica-friendly format

    print output
