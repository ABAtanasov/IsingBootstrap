# -----------------------------------------------------------------
# mixed_ising.py
#
# (Adopted from T. Ohtsuki's cboot example of the same name)
#
# This module uses sage to build the conformal block tables to be
# fed into sdpb, as well as all other relevant sdpb parameters
#
# Depending on the params fed in, this either reads from an in_file
# or constructs the set of points to loop over by using
# point_generator.py
#
# Further depending on the params, this loop will either be run
# directly on the points fed in, or a more complicated method such
# as evelope or bisection will be run
# See 'envelope_method.py' and 'bisection_method.py' respectively
# -----------------------------------------------------------------

import sage.cboot as cb
from sage.misc.cachefunc import cached_function
from subprocess import Popen, PIPE
import numpy as np
import sys
import parser_tools
import envelope_method as env
import bisection_method as bis
from point_generator import generate_points, generate_from_file
from printing_tools import *

mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

sdpb = "./sdpb"
sdpbparams = ["--findPrimalFeasible",
              "--maxRuntime=150000",
              "--findDualFeasible",
              "--noFinalCheckpoint",
              "--dualErrorThreshold","1e-10",
              "--primalErrorThreshold","1e-120"]


context = None
lmax = None
nu_max = None
name = None
keepxml = None
print_sdpb = None
profile = None
odd_scalar_gap = 3
even_scalar_gap = 3
spin_2_gap = None


@cached_function
def prepare_g_0(spin, Delta=None):
    return context.approx_cb(nu_max, spin)


# g_\Delta l
@cached_function
def prepare_g_se(spin, Delta_se, Delta=None):
    g_se = context.approx_cb(nu_max, spin, Delta_1_2=Delta_se, Delta_3_4=Delta_se)
    return g_se


@cached_function
def prepare_g_es(spin, Delta_se, Delta=None):
    g_es = context.approx_cb(nu_max, spin, Delta_1_2=-Delta_se, Delta_3_4=Delta_se)
    return g_es


def prepare_g(spin, Delta_se, Delta=None):
    if Delta == None:
        return (prepare_g_0(spin),
                prepare_g_se(spin, Delta_se),
                prepare_g_es(spin, Delta_se))
    else:
        g_0 = context.gBlock(spin, Delta, 0, 0)
        if not (Delta == 0 and spin == 0):
            g_se = context.gBlock(spin, Delta, Delta_se, Delta_se)
            g_es = context.gBlock(spin, Delta, -Delta_se, Delta_se)
        else:
            g_se = None
            g_es = None
        return (g_0, g_se, g_es)


def make_F(deltas, sector, spin, gap_dict, Delta=None):
    delta_s = context(deltas[0])
    delta_e = context(deltas[1])
    Delta_se = delta_s-delta_e
    if Delta == None:
        try:
            shift=context(gap_dict[(sector, spin)])
        except KeyError:
            if spin == 0:
                shift = context.epsilon
            else:
                shift = 2*context.epsilon + spin
        gs = [x.shift(shift) for x in prepare_g(spin, Delta_se, Delta=Delta)]
    else:
        gs = prepare_g(spin, Delta_se, Delta=Delta)

    if sector == "even":
        F_s_s = context.dot(context.F_minus_matrix(delta_s), gs[0])
        F_e_e = context.dot(context.F_minus_matrix(delta_e), gs[0])

        F_s_e = context.dot(context.F_minus_matrix((delta_s+delta_e)/2), gs[0])
        H_s_e = context.dot(context.F_plus_matrix((delta_s+delta_e)/2), gs[0])
        return [[[F_s_s, 0],
                 [0, 0]],
                [[0, 0],
                 [0, F_e_e]],
                [[0, 0],
                 [0, 0]],
                [[0, F_s_e/2],
                 [F_s_e/2, 0]],
                [[0, H_s_e/2],
                 [H_s_e/2, 0]]]

    elif sector == "odd+":
        F_s_e = context.dot(context.F_minus_matrix((delta_s+delta_e)/2), gs[1])
        F_e_s = context.dot(context.F_minus_matrix(delta_s), gs[2])
        H_e_s = context.dot(context.F_plus_matrix(delta_s), gs[2])

        return [0, 0, F_s_e, F_e_s, -H_e_s]

    elif sector == "odd-":
        F_s_e = context.dot(context.F_minus_matrix((delta_s+delta_e)/2), gs[1])
        F_e_s = context.dot(context.F_minus_matrix(delta_s), gs[2])
        H_e_s = context.dot(context.F_plus_matrix(delta_s), gs[2])

        return [0, 0, -F_s_e, F_e_s, -H_e_s]
    else:
        raise RuntimeError("unknown sector name")


def make_SDP(deltas):
    pvms = []
    if spin_2_gap is not None:
        gaps = {("even", 0): even_scalar_gap, ("odd+", 0): odd_scalar_gap, ("even", 2): spin_2_gap}
    else:
        gaps = {("even", 0): even_scalar_gap, ("odd+", 0): odd_scalar_gap}
    for spin in range(0, lmax):
        if not spin % 2:
            pvms.append(make_F(deltas, "even", spin, gaps))
            pvms.append(make_F(deltas, "odd+", spin, gaps))
        else:
            pvms.append(make_F(deltas, "odd-", spin, gaps))

    epsilon_contribution = make_F(deltas, "even", 0, {}, Delta=deltas[1])
    sigma_contribution = make_F(deltas, "odd+", 0, {}, Delta=deltas[0])
    for m, x in zip(epsilon_contribution, sigma_contribution):
        m[0][0] += x
    if len(deltas) == 3:
        theta = deltas[2]
        V = epsilon_contribution
        constraint = []
        sn = np.sin(theta)
        cs = np.cos(theta)
        for m in V:
            ans = m[0][0] * (cs ** 2) + m[0][1] * sn * cs + m[1][0] * sn * cs + m[1][1] * (sn ** 2)
            constraint.append(ans)
        pvms.append(constraint)
    else:
        pvms.append(epsilon_contribution)
    norm = []
    for v in make_F(deltas, "even", 0, {}, Delta=0):
        norm.append(v[0][0]+v[0][1]+v[1][0]+v[1][1])
    obj = 0
    return context.sumrule_to_SDP(norm, obj, pvms)


def check(deltas, f_out=None):
    # Use cboot to construct sdpb input data:
    start = time.time()
    prob = make_SDP(deltas)
    end = time.time()
    cboot_duration = end - start

    # Write sdpb data to a .xml file in scratch:
    xmlfile = os.path.join(scratchpath, name + ".xml")
    prob.write(xmlfile)
    sdpbargs = [sdpb, "-s", xmlfile] + sdpbparams

    # Run sdpb as a subprocess using Popen,
    # saving the sdpb output and error messages into 'out' and 'err'
    start = time.time()
    out, err = Popen(sdpbargs, stdout=PIPE, stderr=PIPE).communicate()
    end = time.time()
    sdpb_duration = end - start

    # By default, we remove the .xml afterwards:
    if not keepxml:
        os.remove(xmlfile)
    # Optionally, we print the whole sdpb output to the out file:
    if print_sdpb:
        print_out(out, f_out=f_out)
    # If an sdpb error is encountered, print it for the user to see in the out file:
    if err:
        print_err("An error occurred: ", err, f_out=f_out)
        os.system("rm scratch/{}.ck".format(name))
        return False

    durations = (cboot_duration, sdpb_duration)

    # Print all the data of the run (depending on various parameters) to the out file:
    excluded = print_point(deltas, out, name, durations, profile, f_out=f_out)
    return excluded

if __name__ == "__main__":
    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("sage {} -h".format(os.path.abspath(__file__)))
        exit(0)

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("python {} -h".format(os.path.abspath(__file__)))
        exit(0)

    # Parse lines from command line
    parser = parser_tools.build_parser()

    # Take the args as dictionary
    args = parser.parse_args().__dict__

    # Params fed into sdpb
    sdpb_params = {'Lambda': 11, 'lmax': 20, 'nu_max': 8, 'precision': 400, 'maxIters': 500, 'threads': 4}

    # Params fed into the cluster/mixed_ising
    job_params = {'name': "untitled",
                  'res': [1, 1], 'theta_res': 1,
                  'range': None, 'theta_range': None,
                  'dist': None, 'theta_dist': None,
                  'origin': None,
                  'keepxml': False, 'print_sdpb': False, 'profile': False, 'envelope': False,
                  'odd_scalar_gap': 3, 'even_scalar_gap': 3, 'spin_2_gap': None,
                  'max_bisections': None, 'side': "exterior",
                  'sig_spacing': None, 'eps_spacing': None, 'theta_spacing': None,
                  'out_file': False, 'in_file': None}

    # params fed into sdpb
    for key in sdpb_params.keys():
        if args[key] is not None:
            sdpb_params[key] = args[key]
        elif key in ['Lambda', 'lmax', 'nu_max']:
            print "Warning, {} not specified. Using {} = {}.".format(key, key, sdpb_params[key])

    # params characterizing the job
    for key in job_params.keys():
        if args[key] is not None:
            job_params[key] = args[key]

    name = job_params['name']
    keepxml = job_params['keepxml']
    print_sdpb = job_params['print_sdpb']
    profile = job_params['profile']
    envelope = job_params['envelope']

    # Decide whether we print to a file or just print to STDOUT
    f_out = None
    if job_params['out_file']:
        f_out = open("out_files/{}.out".format(name), 'w')

    # Get the points to loop over:
    if job_params['in_file'] is not None:
        f_in = open(job_params['in_file'], 'r')
        points = generate_from_file(job_params, f_in=f_in, f_out=f_out)
    else:
        points = generate_points(job_params, f_out=f_out)

    Lambda = sdpb_params['Lambda']
    lmax = sdpb_params['lmax']
    nu_max = sdpb_params['nu_max']
    odd_scalar_gap = job_params['odd_scalar_gap']
    even_scalar_gap = job_params['even_scalar_gap']
    spin_2_gap = job_params['spin_2_gap']


    sdpbparams.append("--precision={}".format(sdpb_params['precision']))
    sdpbparams.append("--maxThreads={}".format(sdpb_params['threads']))
    sdpbparams.append("--maxIterations={}".format(sdpb_params['maxIters']))

    # Make the context for a 3D conformal field theory
    context = cb.context_for_scalar(epsilon=0.5, Lambda=Lambda)

    assert len(points) > 0

    if job_params['max_bisections'] is not None:
        if job_params['theta_spacing'] is not None:
            spacing = [job_params['sig_spacing'], job_params['eps_spacing'], job_params['theta_spacing']]
        else:
            spacing = [job_params['sig_spacing'], job_params['eps_spacing']]
        bis.bisect_loop(check, points,
                        spacing=spacing, max_bisections=job_params['max_bisections'],
                        f=f_out, side=job_params['side'])
    elif len(points[0]) == 2:
        if envelope:
            env.envelope_loop2D(check, points, f=f_out)
        else:
            for point in points:
                check((point[0], point[1]), f_out=f_out)
    elif len(points[0]) == 3:
        if envelope:
            env.envelope_loop3D(check, points, f=f_out)
        else:
            for point in points:
                check((point[0], point[1], point[2]), f_out=f_out)
    else:
        raise NotImplementedError('Inputs of more than than three CFT data are not yet implemented')


    if f_out is not None:
        f_out.close()

    # TODO: make a 'clean' method at the end, to remove possible lurking .ck files
