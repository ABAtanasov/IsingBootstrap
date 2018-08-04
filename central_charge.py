# -----------------------------------------------------------------
# central_charge.py
#
# Aaron's script (based off of T. Ohtsuki's mixed_ising.py)
# for obtaining the bounds on central charge
# -----------------------------------------------------------------

import sage.cboot as cb
from sage.misc.cachefunc import cached_function
from subprocess import Popen, PIPE
import re
import numpy as np
import sys
import os
import argparse
import time

mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

sdpb = "./sdpb"
sdpbparams = [
              "--noFinalCheckpoint",
              "--dualErrorThreshold=1e-30", "--primalErrorThreshold=1e-30", "--maxRuntime=400000", "--maxComplementarity=1e220"]

context = None
lmax = None
nu_max = None
name = None
keepxml = None
print_sdpb = None

#
@cached_function
def prepare_g_0(spin,Delta=None):
    return context.approx_cb(nu_max,spin)


# g_\Delta l
@cached_function
def prepare_g_se(spin,Delta_se,Delta=None):
    g_se=context.approx_cb(nu_max,spin,Delta_1_2=Delta_se,Delta_3_4=Delta_se)
    return g_se


@cached_function
def prepare_g_es(spin,Delta_se,Delta=None):
    g_es=context.approx_cb(nu_max,spin,Delta_1_2=-Delta_se,Delta_3_4=Delta_se)
    return g_es


def prepare_g(spin,Delta_se,Delta=None):
    if Delta==None:
        return (prepare_g_0(spin),
                prepare_g_se(spin,Delta_se),
                prepare_g_es(spin,Delta_se))
    else:
        g_0=context.gBlock(spin,Delta,0,0)
        if not (Delta==0 and spin==0):
            g_se=context.gBlock(spin,Delta,Delta_se,Delta_se)
            g_es=context.gBlock(spin,Delta,-Delta_se,Delta_se)
        else:
            g_se=None
            g_es=None
        return (g_0,g_se,g_es)


def make_F(deltas,sector,spin,gap_dict,Delta=None):
    delta_s=context(deltas[0])
    delta_e=context(deltas[1])
    Delta_se=delta_s-delta_e
    if Delta==None:
        try:
            shift=context(gap_dict[(sector,spin)])
        except KeyError:
            if spin==0:
                shift=context.epsilon
            else:
                shift=2*context.epsilon+spin
        gs=[x.shift(shift) for x in prepare_g(spin,Delta_se,Delta=Delta)]
    else:
        gs=prepare_g(spin,Delta_se,Delta=Delta)

    if sector=="even":
        F_s_s=context.dot(context.F_minus_matrix(delta_s),gs[0])
        F_e_e=context.dot(context.F_minus_matrix(delta_e),gs[0])

        F_s_e=context.dot(context.F_minus_matrix((delta_s+delta_e)/2),gs[0])
        H_s_e=context.dot(context.F_plus_matrix((delta_s+delta_e)/2),gs[0])
        return [[[F_s_s,0],
                   [0,0]],
                  [[0,0],
                  [0,F_e_e]],
                  [[0,0],
                   [0,0]],
                  [[0,F_s_e/2],
                   [F_s_e/2,0]],
                  [[0,H_s_e/2],
                   [H_s_e/2,0]]]

    elif sector=="odd+":
        F_s_e=context.dot(context.F_minus_matrix((delta_s+delta_e)/2),gs[1])
        F_e_s=context.dot(context.F_minus_matrix(delta_s),gs[2])
        H_e_s=context.dot(context.F_plus_matrix(delta_s),gs[2])

        return [0,0,F_s_e,F_e_s,-H_e_s]

    elif sector=="odd-":
        F_s_e=context.dot(context.F_minus_matrix((delta_s+delta_e)/2),gs[1])
        F_e_s=context.dot(context.F_minus_matrix(delta_s),gs[2])
        H_e_s=context.dot(context.F_plus_matrix(delta_s),gs[2])

        return [0,0,-F_s_e,F_e_s,-H_e_s]
    else: raise RuntimeError("unknown sector name")


def make_SDP(deltas, theta=None):

    pvms=[]
    gaps={("even",0):3,("odd+",0):3, ("even", 2):3}
    for spin in range(0,lmax):
        if not spin%2:
            pvms.append(make_F(deltas,"even",spin,gaps))
            pvms.append(make_F(deltas,"odd+",spin,gaps))
        else:
            pvms.append(make_F(deltas,"odd-",spin,gaps))

    epsilon_contribution=make_F(deltas,"even",0,{},Delta=deltas[1])
    sigma_contribution=make_F(deltas,"odd+",0,{},Delta=deltas[0])
    for m,x in zip(epsilon_contribution,sigma_contribution):
        m[0][0]+=x
    if theta is not None:
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
    stress_tensor = make_F(deltas, "even", 2, {}, Delta = 3)
    pvms.append(stress_tensor)
    norm=[]
    for v in make_F(deltas,"even",0,{},Delta=0):
        norm.append(v[0][0]+v[0][1]+v[1][0]+v[1][1])
    obj=0
    return context.sumrule_to_SDP(norm, obj, pvms)

#we use this function if we're minimizing c instead of using trivial objective function
def make_SDP_cmin(deltas, theta=None):
    delta_s = context(deltas[0])
    delta_e = context(deltas[1])
    pvms=[]
    gaps={("even",0):epspgap,("odd+",0):sigpgap}
    for spin in range(0,lmax):
        if not spin%2:
            pvms.append(make_F(deltas,"even",spin,gaps))
            pvms.append(make_F(deltas,"odd+",spin,gaps))
        else:
            pvms.append(make_F(deltas,"odd-",spin,gaps))

    epsilon_contribution=make_F(deltas,"even",0,{},Delta=deltas[1])
    sigma_contribution=make_F(deltas,"odd+",0,{},Delta=deltas[0])
    for m,x in zip(epsilon_contribution,sigma_contribution):
        m[0][0]+=x
    if theta is not None:
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
    norm=[]
    obj=[]
    stress_tensor = make_F(deltas, "even", 2, {}, Delta = 3)
    rat00 = (delta_s/delta_e)
    rat11 = (delta_e/delta_s)
    for v in stress_tensor:
        ans = v[0][1] + v[1][0] + v[0][0] * rat00 + v[1][1] * rat11
        norm.append(ans)
    for v in make_F(deltas,"even",0,{},Delta=0):
        obj.append(v[0][0]+v[0][1]+v[1][0]+v[1][1])
    return context.sumrule_to_SDP(norm, obj, pvms)


def check(deltas, checkpoint, theta=None, f=None):


    # if not checkpointing, also make the problem for SDPB
    checkpointpath = os.path.join(scratchpath, name + ".xml")
    if not os.path.exists(checkpointpath):
        f.write("SDPB...\n")
        f.flush()
        time.sleep(4)
        prob = make_SDP_cmin(deltas, theta)

        xmlfile = os.path.join(scratchpath, name + ".xml")
        ckfile = os.path.join(scratchpath, name + ".ck")
        prob.write(xmlfile)
        sdpbargs = [sdpb, "-s", xmlfile] + sdpbparams
    else:
        xmlfile = os.path.join(scratchpath, name + ".xml")
        sdpbargs = [sdpb, "-s", xmlfile] + sdpbparams



    out, err = Popen(sdpbargs, stdout=PIPE, stderr=PIPE).communicate()
    sol = re.compile(r'primalObjective *= *([^ ]+) *$', re.MULTILINE) \
        .search(out).groups()[0]
    f.write("Central Charge Bound for ({},{}, {}): {}\n".format(deltas[0], deltas[1], theta,
                                                                -(deltas[0] * deltas[1]) / float(sol)))






if __name__ == "__main__":

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("sage {} -h".format(os.path.abspath(__file__)))
        exit(0)

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("python {} -h".format(os.path.abspath(__file__)))
        exit(0)

    parser = argparse.ArgumentParser()

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

    # --------------------------------------
    # Args for sdpb
    # --------------------------------------
    parser.add_argument("-p", "--precision", type=int,
                        help="working precision for sdpb calculations")
    parser.add_argument("-i", "--maxIters", type=int,
                        help="max number of sdpb iterations")
    parser.add_argument("--threads", type=int,
                        help="maximum threads used by OpenMP")

    #point testing
    parser.add_argument("-ds", "--delta_s", type=float,
                        help="delta sigma")
    parser.add_argument("-de", "--delta_e", type=float,
                        help="delta epsilon")
    parser.add_argument("-th", "--theta", type=float)
    parser.add_argument("-cp", "--checkpoint", type=str, help="are we starting from an SDPB checkpoint?")
    parser.add_argument("-spg", "--sigpgap", type=float,
                        help="delta sigma")
    parser.add_argument("-epg", "--epspgap", type=float,
                        help="delta epsilon")

    args = parser.parse_args().__dict__

    args['checkpoint'] = (args['checkpoint'] == 'True')

    # Params fed into sdpb
    sdpb_params = {'Lambda': 11, 'lmax': 20, 'nu_max': 8, 'precision': 400, 'maxIters': 500, 'threads': 4, 'checkpoint':False}

    # Params fed into the cluster/mixed_ising
    job_params = {'name':"untitled",
            'res':[1, 1], 'theta_res':1,
            'range': None, 'theta_range': None,
            'dist': None, 'theta_dist':None,
            'origin': None,
            'keepxml':False, 'print_sdpb':False,
            'out_file':False, 'in_file': None, 'delta_s':0, 'delta_e':0, 'theta':None, 'sigpgap':3.0, 'epspgap': 3.0}

    # params fed into sdpb
    for key in sdpb_params.keys():
        if args[key]:
            sdpb_params[key] = args[key]
        elif key in ['Lambda', 'lmax', 'nu_max']:
            print "Warning, {} not specified. Using {} = {}.".format(\
                    key, key, sdpb_params[key])


    # params characterizing the job
    for key in job_params.keys():
        if args[key]:
            job_params[key]=args[key]

    name = job_params['name']
    keepxml = job_params['keepxml']
    print_sdpb = job_params['print_sdpb']
    delta_s = job_params['delta_s']
    delta_e = job_params['delta_e']
    theta = job_params['theta']
    sigpgap = job_params['sigpgap']
    epspgap = job_params['epspgap']
    # Decide whether we print to a file or just print out
    f_out = None
    if job_params['out_file']:
        f_out = open("out_files/{}.out".format(name), "a+")



    Lambda = sdpb_params['Lambda']
    lmax = sdpb_params['lmax']
    nu_max = sdpb_params['nu_max']
    checkpoint = sdpb_params['checkpoint']

    sdpbparams.append("--precision={}".format(sdpb_params['precision']))
    sdpbparams.append("--maxThreads={}".format(sdpb_params['threads']))
    sdpbparams.append("--maxIterations={}".format(sdpb_params['maxIters']))

    context=cb.context_for_scalar(epsilon=0.5,Lambda=Lambda)

    check((delta_s, delta_e), checkpoint, theta, f=f_out)


    f_out.close()
    # make a 'clean' method at the end, to remove possibly lurking .ck files
