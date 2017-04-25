# -----------------------------------------------------------------
# point_generator.py
# 
# -----------------------------------------------------------------

import numpy as np
import os
import re

mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

# --------------------------------------------------------
# Makes a numpy array from a to b with a given resolution
# --------------------------------------------------------
def mkrange(a, b, resolution):
    if resolution == 1:
        return np.array([0.5*(a+b)])
    else:
        return np.linspace(a, b, num = resolution)


# --------------------------------------------------------
# Prints to the file f or stdout if no f is specified
# --------------------------------------------------------
def print_out(string, f=None):
    if f is not None:
        f.write(string)
        f.write("\n")
    else:
        print string

# --------------------------------------------------------
# If f_in is specified, generates the points from that
# in_file
#
# Otherwise, uses params to construct the correct loop
# --------------------------------------------------------
def generate_points(params, f_in=None, f_out=None):

    if f_in is not None:
        return generate_from_file(params, f_in=f_in, f_out=f_out)

    Dsig = 0.518154
    Deps = 1.41267
    distance = (0.002, 0.02)
    theta0 = 0.969260330903202

    if params['origin'] and not params['range']:
        Dsig, Deps = params['origin']
    if params['dist'] and not params['range']:
        dist = float(params['dist'])
        distance = (dist, 10 * dist)
    sig_min = Dsig - distance[0]
    sig_max = Dsig + distance[0]
    eps_min = Deps - distance[1]
    eps_max = Deps + distance[1]
    if params['range']:
        sig_min, sig_max, eps_min, eps_max = params['range']
    sig_res, eps_res = params['res']

    use_theta = False
    if params['theta_range']:
        use_theta = True
        theta_min, theta_max = params['theta_range']
    elif params['theta_dist']:
        use_theta = True
        theta_min = theta0 - params['theta_dist']
        theta_max = theta0 + params['theta_dist']
    theta_res = params['theta_res']

    print_out("Using {}".format(params), f=f_out)
    print_out("with resolutions = ({}, {}), ".format(sig_res, eps_res) \
          + "Delta window = (({}, {}), ({}, {})), ".format(
        sig_min, sig_max, eps_min, eps_max), f=f_out)

    if use_theta:
        print_out("use_theta window = ({}, {}), with resolution {}".format(
            theta_min, theta_max, theta_res), f=f_out)

    sigmas = mkrange(sig_min, sig_max, sig_res)
    epsilons = mkrange(eps_min, eps_max, eps_res)
    points = []
    if use_theta:
        thetas = mkrange(theta_min, theta_max, theta_res)

    for sigma in sigmas:
        for epsilon in epsilons:
            if use_theta:
                for theta in thetas:
                    points.append([sigma, epsilon, theta])
            else:
                points.append([sigma, epsilon])

    return points

# --------------------------------------------------------
# Given an in_file, prints the params to an out_file
# (or stdout by default), and returns the list of points
# contained in that file
# --------------------------------------------------------
def generate_from_file(params, f_in, f_out):

    print_out("Using {}".format(params), f=f_out)
    print_out("from file {}".format(f_in.name), f=f_out)

    number_data = re.compile("[\d]+.[\d]*")
    points = []

    # The following gives the list of points contained in the file
    for line in f_in:
        data = number_data.findall(line)
        if data:
            points.append(map(lambda x: float(x), data))

    return points

# --------------------------------------------------------
# Generates a set of points, either from params or an
# in_file, to and prints these to an out_file
# corresponding to the name of the job
#
# If multiple batches are specified, this divides the
# points into equal batches to be written to seperate
# files, one for each job that will be called, and
# labelled in the form:
#
#        <name>_<batch #>of<total batches>.pts
#
# --------------------------------------------------------
def generate_to_file(params, batches=1, f_in=None):

    points = generate_points(params, f_in=f_in)
    name = params['name']
    num_points = len(points)

    for batch in range(batches):
        f_new = open("scratch/{}_{}of{}.pts".format(name, batch + 1, batches), 'w')

        begin = (batch * num_points)/batches
        end   = ((batch + 1) * num_points)/batches
        for point in points[begin:end]:
            f_new.write("{}\n".format(point))

        f_new.close()

