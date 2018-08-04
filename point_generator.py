# --------------------------------------------------------
# point_generator.py
#
# This script generates and manipulates the arrays of
# points to be looped over using various methods
# (e.g. enveleope, bisection) and fed into sdpb
# --------------------------------------------------------

import numpy as np
from printing_tools import *

mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

# --------------------------------------------------------
# Makes a numpy array from a to b with a given resolution
# --------------------------------------------------------
def mkrange(a, b, resolution):
    if resolution == 1:
        return np.array([0.5*(a+b)])
    else:
        return np.linspace(a, b, num=resolution)


# --------------------------------------------------------
# Converts an array of 3D point values, interpretted as
# tuples: (\Delta_\sigma, \Delta_\epsilon, \theta)
# into a dict of the base tuples
# (\Delta_\sigma, \Delta_\epsilon) as keys
# and the theta values as vals
# --------------------------------------------------------
def array2dict3D(points):
    assert len(points) > 0 and len(points[0]) == 3
    base_points = dict()
    # Build the dictionary:
    for point in points:
        base_point = (point[0], point[1])
        theta = point[2]
        entry = base_points.get(base_point)
        if entry is None:
            base_points[base_point] = [theta]
        else:
            base_points[base_point].append(theta)

    # Sort the dictionary
    for thetas in base_points.values():
        thetas.sort()

    return base_points

# --------------------------------------------------------
# Converts an array of 2D point values, interpretted as
# tuples: (\Delta_\sigma, \Delta_\epsilon)
# into a dict of the sigma values as keys
# and the epsilon values as vals
# --------------------------------------------------------
def array2dict2D(points):
    assert len(points) > 0 and len(points[0]) == 2
    base_points = dict()
    # Build the dictionary:
    for point in points:
        base_point = (point[0], )
        eps = point[1]
        entry = base_points.get(base_point)
        if entry is None:
            base_points[base_point] = [eps]
        else:
            base_points[base_point].append(eps)

    # Sort the dictionary
    for epsilons in base_points.values():
        epsilons.sort()

    return base_points

# --------------------------------------------------------
# Given an in_file, prints the params to an out_file
# (or STDOUT by default), and returns the list of points
# contained in that file
# --------------------------------------------------------
def generate_from_file(params=None, f_in=None, f_out=None):

    if f_out is not None and f_in is not None:
        print_out("Using {}".format(params), f_out=f_out)
        print_out("from file {}".format(f_in.name), f_out=f_out)

    number_data = re.compile("-?[\d]+.[\d]*")
    points = []

    # The following gives the list of points contained in the file
    for line in f_in:
        data = number_data.findall(line)
        if data:
            points.append(map(lambda x: float(x), data))

    return points


# --------------------------------------------------------
# If f_in is specified, generates the points from that
# in_file
#
# Otherwise, uses params to construct the correct array
# of points
# --------------------------------------------------------
def generate_points(params, f_in=None, f_out=None):

    # Usually we will just read in from a file (can effectively ignore things below this)
    if f_in is not None:
        return generate_from_file(params, f_in=f_in, f_out=f_out)

    # Otherwise I will assume we are looking to resolve information near the 3D Ising CFT point
    Dsig = 0.518154
    Deps = 1.41267
    distance = (0.002, 0.02)
    theta0 = 0.969260330903202

    # Or some other origin/distance specified by the user
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

    print_out("Using {}".format(params), f_out=f_out)
    print_out("with resolutions = ({}, {}), ".format(sig_res, eps_res) \
              + "Delta window = (({}, {}), ({}, {})), ".format(
        sig_min, sig_max, eps_min, eps_max), f_out=f_out)

    if use_theta:
        print_out("use_theta window = ({}, {}), with resolution {}".format(
            theta_min, theta_max, theta_res), f_out=f_out)

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

    write_update(f_out)
    return points


# --------------------------------------------------------
# Generates a set of points, either from params or an
# in_file, to and prints these to an out_file
# corresponding to the name of the job
#
# If multiple batches are specified, this divides the
# points into equal batches to be written to separate
# files, one for each job that will be called, and
# labelled in the form:
#
#        <name>_<batch #>of<total batches>.pts
#
# --------------------------------------------------------
def generate_to_file(params, batches=1, f_in=None):

    points = generate_points(params, f_in=f_in)
    name = params['name']

    assert len(points) > 0
    assert len(points[0]) == 2 or len(points[0]) == 3

    # In the case of an envelope method this bundles the submit files by
    # Base point in the plane, so that the full fiber is accessible for a given job
    if params['envelope']:
        if len(points[0]) == 2:
            base_points = array2dict2D(points)
        elif len(points[0]) == 3:
            base_points = array2dict3D(points)
        num_points = len(base_points.keys())
        batch = 0
        point_num = 0
        f_new = open("scratch/{}_{}of{}.pts".format(name, batch + 1, batches), 'w')
        for base_point, fiber in base_points.iteritems():
            if point_num > num_points/batches and batch+1 != batches:
                batch += 1
                point_num = 1
                f_new.close()
                f_new = open("scratch/{}_{}of{}.pts".format(name, batch + 1, batches), 'w')
            for value in fiber:
                if len(points[0]) == 2:
                    f_new.write("({}, {})\n".format(base_point[0], value))
                elif len(points[0]) == 3:
                    f_new.write("({}, {}, {})\n".format(base_point[0], base_point[1], value))
            point_num += 1

    # In the case of a 2D scan over scaling dimensions, we break down
    # the batches equally.
    else:
        num_points = len(points)
        for batch in range(batches):
            f_new = open("scratch/{}_{}of{}.pts".format(name, batch + 1, batches), 'w')

            begin = (batch * num_points)/batches
            end   = ((batch + 1) * num_points)/batches
            for point in points[begin:end]:
                f_new.write("{}\n".format(point))


    f_new.close()

