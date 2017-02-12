import os
import math
import sys
import numpy as np
import argparse

def mkrange(a, b, resolution):
    if resolution == 1:
        return np.array([0.5*(a+b)])
    else:
        return np.linspace(a, b, num = resolution)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-B", "--batches", type=str, nargs='+',\
            help="info for how jobs are submitted into batches")

    parser.add_argument("-N", "--name", type = str,\
            help="name for the associated files")
    parser.add_argument("-P", "--program", type = str,\
            help="name of the program to execute")
    parser.add_argument("-L","--Lambda", type = int, \
           help="maximum derivative order")
    parser.add_argument("-l", "--lmax", type = int, \
            help="angular momentum cutoff")
    parser.add_argument("-nu", "--nu_max", type = int, \
            help="maximum number of poles")
    parser.add_argument("-p", "--precision", type =int, \
            help="working precision for calculations")
    parser.add_argument("--res", type = int, nargs = 2,\
            help="number of sampling points along each axis")
    parser.add_argument("--theta_res", type = int, nargs = 1,\
            help="number of sampling points over theta")
    parser.add_argument("--range", type = float, nargs = 4,\
            help="4 floats xmin xmax ymin ymax")
    parser.add_argument("--theta_range", type = float, nargs = 2,\
            help="2 floats theta_min theta_max")
    parser.add_argument("--mem", type = int,\
            help="maximum memory allocated per node in cluster")
    parser.add_argument("--ndays", type = int,\
            help="number of days to run process on cluster")
    parser.add_argument("--threads", type = int, \
            help="maximum threads used by OpenMP")
    parser.add_argument("-q", "--queue", type = str,\
            help="queue to submit to")
    args = parser.parse_args().__dict__

    if len(sys.argv) == 1:
        os.system("python {} -h".format(os.path.abspath(__file__)))
        exit(0)

    if args["program"]==None:
        print "Error: program must be specified"
        exit(0)

    if args["batches"] == None:
        print "Error: no batch division specified"
        exit(0)

    else: batches = args["batches"]

    cmd = "python submit.py"

    modified_keys = ["name", "batches", "res", "theta_res", "range", "theta_range"]

    if len(batches) < 1 or len(batches) > 2:
        print "Error: --batch flag should have 1 or 2 arguments."


    res = [1, 1]
    theta_res = 1
    range = [0.518154, 0.518154, 1.41267, 1.41267]
    theta_range = [0.96926, 0.96926]

    if args["res"]: res = args["res"]
    if args["theta_res"]: theta_res = args["theta_res"]
    if args["range"]: range = args["range"]
    if args["theta_range"]: theta_range = args["theta_range"]

    for key in args:
        if not (key in modified_keys) :
            cmd += " --" + key + " " + args[key]

    if len(batches) == 1:
        theta_bits = mkrange(theta_range[0], theta_range[1], batches+1)
        theta_ranges = [[theta_bits[i+1], theta_bits[i]] for i in range(batches)]
        for batch in range(batches):
            additional_info = " --res {} {}".format(res[0], res[1]) \
                    +" --theta_res {}".format(theta_res/batches)\
                    +" --range {} {} {} {}".format(\
                    range[0], range[1], range[2], range[3])\
                    +" --theta_range {} {}".format(\
                    theta_ranges[batch][0], theta_ranges[batch][1])
            os.system(cmd + additional_info)

    elif len(batches) == 2:
        sig_bits = mkrange(range[0], range[1], batches[0]+1)
        eps_bits = mkrange(range[2], range[3], batches[1]+1)
        sig_ranges = [[sig_bits[i+1], sig_bits[i]] for i in range(batches[0])]
        eps_ranges = [[eps_bits[i+1], eps_bits[i]] for i in range(batches[1])]
        for sig_batch in range(batches[0]):
            for eps_batch in range(batches[1]):
                additional_info = " --res {} {}".format(\
                        res[0]/batches[0], res[1]/batches[1])\
                        + " --range {} {} {} {}".format(\
                        sig_ranges[sig_batch][0], sig_ranges[sig_batch][1],\
                        eps_ranges[eps_batch][0], eps_ranges[eps_batch][1])
                if args["program"] == "theta_scan":
                    additional_info += " --theta_res {}".format(theta_res)\
                            + " --theta_range {} {}".format(\
                            theta_range[0], theta_range[1])
                os.system(cmd + additional_info)





