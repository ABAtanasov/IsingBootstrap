# -----------------------------------------------------------------
# spectrum_data_manipulator.py
#
# Given a mathematica file resulting from running spectrum.py,
# Parses it into the format we use for reading out spectral data
# -----------------------------------------------------------------

import os
import argparse
import subprocess

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--name", type=str)
    parser.add_argument("--Tgap", type=float)
    parser.add_argument("--oddgap", type=float)
    parser.add_argument("--evengap", type=float)
    parser.add_argument("--maxOperators", type=int)

    args = parser.parse_args().__dict__

    #gap information so we know what to add spectrum extraction data to if not the unitarity bound
    name = args['name']
    Tgap = args['Tgap']
    oddgap = args['oddgap']
    evengap = args['evengap']
    maxOperators = args['maxOperators']
    gaps = [evengap, oddgap, 2, Tgap, 3, 4, 5, 5, 6, 7, 7]

    #getting the central charge data in python array format
    spectrumdata = eval(subprocess.check_output(["python","ccout.py", name]))
    print "Length of spectrum data: ", len(spectrumdata)
    for i in range(0, len(spectrumdata)):
        # Open the file path for the spectrum data obtained from running spectrum.py
        spectrumpath = "scratch/{}_{}.spectrum.m".format(name, i)
        if os.path.exists(spectrumpath):
            spectrumfile = open(spectrumpath, "r")
            lines = spectrumfile.readlines()
            operatordata = []
            for line in lines:
                startindex = line.index('>')+2
                endindex = len(line)-2
                trimmed = line[startindex:endindex]
                trimmed = str(trimmed)
                trimmed = trimmed.replace("{", "[")
                trimmed = trimmed.replace("}", "]")
                data = eval(trimmed)
                operatordata.append(data)
            for j in range(0, min(maxOperators, len(operatordata))):
                if len(operatordata[j]) == 1:
                    spectrumdata[i].append("null")
                if len(operatordata[j]) > 0:
                    dims = [operatordata[j][0][0]+gaps[j], operatordata[j][1][0]+gaps[j], operatordata[j][2][0]+gaps[j]]
                    spectrumdata[i].append(dims)
                elif len(operatordata[j]) > 1:
                    shift = operatordata[j][1][0]
                    spectrumdata[i].append(shift+gaps[j])


    spectrumdata = str(spectrumdata)
    spectrumdata = spectrumdata.replace("[", "{")
    spectrumdata = spectrumdata.replace("]", "}")
    print spectrumdata
    #now we append in order the Delta sigma', Delta epsilon'

