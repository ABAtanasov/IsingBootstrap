# --------------------------------------------------------
# parser_tools.py
#
# All of the argument parser setup is handled by the
# method 'build_parser' of this module
# --------------------------------------------------------

import argparse


def build_parser():
    parser = argparse.ArgumentParser()

    # --------------------------------------
    # Args for submit.py
    # --------------------------------------
    parser.add_argument("-B", "--batches", type=int,
                        help="info for how jobs are submitted into batches")
    parser.add_argument("-f", "--file", type=str,
                        help="optional filename from which we read the points")

    # --------------------------------------
    # Args for mixed_ising.py
    # --------------------------------------
    parser.add_argument("-N", "--name", type=str,
                        help="name for the associated files")
    parser.add_argument("-L", "--Lambda", type=int,
                        help="maximum derivative order")
    parser.add_argument("-l", "--lmax", type=int,
                        help="angular momentum cutoff")
    parser.add_argument("-nu", "--nu_max", type=int,
                        help="maximum number of poles")
    parser.add_argument("--res", type=int, nargs=2,
                        help="number of sampling points along each axis")
    parser.add_argument("--theta_res", type=int,
                        help="number of sampling points over use_theta")
    parser.add_argument("--dist", type=float,
                        help="distance of Delta_sigma window from the 3D Ising point")
    parser.add_argument("--theta_dist", type=float,
                        help="distance of use_theta window from the 3D Ising use_theta")
    parser.add_argument("--range", type=float, nargs=4,
                        help="4 floats xmin xmax ymin ymax")
    parser.add_argument("--origin", type=float, nargs=2,
                        help="2 floats x_origin y_origin")
    parser.add_argument("--theta_range", type=float, nargs=2,
                        help="2 floats theta_min theta_max")
    parser.add_argument("--in_file", type=str,
                        help="file to read points from")
    parser.add_argument("--out_file", type=bool,
                        help="do we print out to a file?")
    parser.add_argument("--keepxml",
                        help="Do we keep the xml? Default is no.")
    parser.add_argument("--print_sdpb",
                        help="Do we print the sdpb output? Default is no.")
    parser.add_argument("--profile",
                        help="Do we profile the time taken? Default is no.")
    parser.add_argument("--envelope",
                        help="Option to use the \'envelope\' method of attack for theta scan")

    # --------------------------------------
    # Args for gap assumptions
    # --------------------------------------
    parser.add_argument("--odd_scalar_gap", type=float,
                        help="Option to change gap assumptions on sigma\'")
    parser.add_argument("--even_scalar_gap", type=float,
                        help="Option to change gap assumptions on Z2 even epsilon\'")
    parser.add_argument("--spin_2_gap", type=float,
                        help="Option to change gap assumptions on Z2 even T\'")

    # --------------------------------------
    # Args for Bisection
    # --------------------------------------
    parser.add_argument("--max_bisections", type=int,
                        help="Maximum number of bisections we run")
    parser.add_argument("--side", type=str,
                        help="Are we starting from the inside or outside of the region during bisection?")
    parser.add_argument("-ssp", "--sig_spacing", type=float,
                        help="initial bisection spacing")
    parser.add_argument("-esp", "--eps_spacing", type=float,
                        help="initial bisection spacing")
    parser.add_argument("-tsp", "--theta_spacing", type=float,
                        help="initial bisection spacing")
    parser.add_argument("-b", "--num_bisections", type=int)

    # --------------------------------------
    # Args for sdpb
    # --------------------------------------
    parser.add_argument("-p", "--precision", type=int,
                        help="working precision for sdpb calculations")
    parser.add_argument("-i", "--maxIters", type=int,
                        help="max number of sdpb iterations")
    parser.add_argument("--threads", type=int,
                        help="maximum threads used by OpenMP")

    # --------------------------------------
    # Args for the cluster
    # --------------------------------------
    parser.add_argument("--mem", type=int,
                        help="maximum memory in GB allocated per node in cluster")
    parser.add_argument("--ndays", type=int,
                        help="number of days to run process on cluster")
    parser.add_argument("-q", "--queue", type=str,
                        help="queue to submit to")

    return parser
