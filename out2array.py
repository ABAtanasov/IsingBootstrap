# -----------------------------------------------------------------
# out2array.py
#
# This module allows us to format the files output from mixed_ising
# into a style that is easy to feed in to mathematica directly.
# -----------------------------------------------------------------

import sys
import re
import os
import numpy as np
from point_generator import generate_from_file

def read_lines(lines, points, form="is not excluded"):
    inclusion = re.compile("(\([-?\d.\d, ]+\)) " + form)
    number_data = re.compile("-?[\d]+.[\d]+")     # general format for a double outputted
    for line in lines:
        inclusionstring = inclusion.search(line)
        if inclusionstring is not None:                              # check that this point is included
            data = number_data.findall(inclusionstring.groups()[0])  # get the 2 or 3 coordinates of the point
            points.append(map(lambda x: float(x), data)) # add it as an array of floats


def read_bisections(lines, points):
    inclusion = re.compile("(\([-?\d.\d, ]+\)) is")  # regex compile for mixed_ising output
    number_data = re.compile("[\d]+.[\d]+")
    line = len(lines)-1
    boundary_found = False
    prev_point = None; prev_point2 = None; point = None
    while line >=0:
        inclusionstring = inclusion.search(lines[line])
        if inclusionstring is None:
            line -= 1   # In case we have errors in certain outputs, we won't break, we'll just keep going
            continue

        data = number_data.findall(inclusionstring.groups()[0]) # find the 2 or 3 coordinates of this point
        prev_point2 = prev_point
        prev_point = point
        point = np.array(map(lambda x: float(x), data))
        if boundary_found and prev_point2 is not None:
            if np.linalg.norm(prev_point - prev_point2) > np.linalg.norm(point - prev_point):
                boundary_found = False
        if not boundary_found and 'not excluded' in lines[line]: # and 'not' not in lines[line]:
            points.append(map(lambda x: round(x, 10), list(point)))  # add it as an array of floats
            boundary_found = True
        line -= 1



if __name__ == "__main__":

    points = []
    envelope = False
    bisect = False
    if len(sys.argv) > 1 and sys.argv[1] == "envelope":
        envelope = True
    if len(sys.argv) > 1 and sys.argv[1] == "bisect":
        bisect = True
    # If one argument is supplied, it'll search through the out files
    # for all those matching that argument
    if len(sys.argv) > 1:
        form = "is not excluded"
        if envelope:
            form = "is excluded"
        for filename in os.listdir("out_files"):
            if (sys.argv[-1]+'_') in filename:
                f_out = open("out_files/{}".format(filename), 'r')
                if bisect:
                    read_bisections(f_out.readlines(), points)
                else:
                    read_lines(f_out.readlines(), points, form=form)
                f_out.close()

    # Otherwise we read from the file that is fed into the program
    else:
        if bisect:
            read_bisections(sys.stdin.readlines(), points)
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
