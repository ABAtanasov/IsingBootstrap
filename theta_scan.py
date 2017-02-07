import sage.cboot as cb
from sage.misc.cachefunc import cached_function
from subprocess import Popen, PIPE
import re
import numpy as np
import sys
import os
import argparse


mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

sdpb="./sdpb"
sdpbparams=["--findPrimalFeasible","--findDualFeasible","--noFinalCheckpoint"]

context = None
lmax   = None
nu_max = None
name = None

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

def make_SDP(deltas, theta):
    pvms=[]
    gaps={("even",0):3,("odd+",0):3}
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

    V = epsilon_contribution
    constraint = []
    sn = np.sin(theta)
    cs = np.cos(theta)
    for m in V:
        ans = m[0][0]*(cs**2) + m[0][1]*sn*cs + m[1][0]*sn*cs + m[1][1]*(sn**2)
        constraint.append(ans)
    pvms.append(constraint)
    norm=[]
    for v in make_F(deltas,"even",0,{},Delta=0):
       norm.append(v[0][0]+v[0][1]+v[1][0]+v[1][1])
    obj=0
    return context.sumrule_to_SDP(norm,obj,pvms)

def check(deltas, theta):
    prob=make_SDP(deltas, theta)

    xmlfile = os.path.join(scratchpath, name+".xml")

    prob.write(xmlfile)
    sdpbargs=[sdpb,"-s",xmlfile]+sdpbparams
    out, err=Popen(sdpbargs,stdout=PIPE,stderr=PIPE).communicate()

    sol = re.compile(r'found ([^ ]+) feasible').search(out)
    sol = sol.groups()[0]
    if sol=="dual":
        print("({}, {}, {}) is excluded."\
            .format(deltas[0], deltas[1], theta))
    elif sol=="primal":
        print("({}, {}, {}) is not excluded."\
            .format(deltas[0], deltas[1], theta))
    else:
        raise RuntimeError

def mkrange(a,b, resolution):
    if resolution == 1:
        return np.array([0.5*(a+b)])
    else:
        return np.linspace(a, b, num = resolution)

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--name", type = str,\
            help="name for the associated files")
    parser.add_argument("-L","--Lambda", type = int, \
           help="maximum derivative order")
    parser.add_argument("-l", "--lmax", type = int, \
            help="angular momentum cutoff")
    parser.add_argument("-nu", "--nu_max", type = int, \
            help="maximum number of poles")
    parser.add_argument("-p", "--precision", type =int, \
            help="working precision for calculations")
    parser.add_argument("--res", type = int,\
            help="number of sampling points along each axis")
    parser.add_argument("--theta_res", type = int,\
            help="number of sampling points over theta")
    parser.add_argument("--dist", type = float,\
            help="distance of Delta_sigma window from the 3D Ising point")
    parser.add_argument("--theta_dist", type = float, \
            help="distance of theta window from the 3D ising theta")
    parser.add_argument("--threads", type = int, \
            help="maximum threads used by OpenMP")
    args = parser.parse_args()

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("sage {} -h".format(os.path.abspath(__file__)))
        exit(0)

    if not args.Lambda:
        print "No Lambda specified."
        exit(1)
    if not args.lmax:
        print "No lmax specified."
        exit(1)
    if not args.nu_max:
        print "No nu_max specified."
        exit(1)

    name   = args.name
    Lambda = args.Lambda
    lmax   = args.lmax
    nu_max = args.nu_max

    res       = 1
    theta_res = 1


    if args.res:
        res = args.res
    if args.theta_res:
        theta_res = args.theta_res

    distance   = (0.002, 0.02)
    theta_dist = 0.1
    if args.dist:
        dist = float(args.dist)
        distance = (dist, 10*dist)
    if args.theta_dist
        theta_dist = float(args.theta_dist)

    precision = 400
    threads = 4
    if args.precision:
        precision = args.precision
        sdpbparams.append("--precision={}".format(args.precision))
    if args.threads:
        threads = args.threads
        sdpbparams.append("--maxThreads={}".format(args.threads))

    print "Using Lambda = {}, lmax = {}, nu_max = {}, precision = {}".format(\
            Lambda, lmax, nu_max, precision)
    print "with resolutions = ({}, {}), ".format(res, theta_res)\
            + "Delta window = ({}, {}), ".format(distance[0], distance[1])\
            + "Theta window = {}, ".format(theta_dist)
            + "threads = {}".format(threads)

    context=cb.context_for_scalar(epsilon=0.5,Lambda=Lambda)

    Dsig = 0.518154
    Deps = 1.41267
    theta0 = 0.969260330903202

    for delta_s in mkrange(Dsig - distance[0], Dsig + distance[0], res):
        for delta_e in mkrange(Deps - distance[1], Deps + distance[1], res):
            for theta in mkrange(theta0 - theta_dist, theta0 + theta_dist, theta_res):
                check((delta_s, delta_e), theta)
