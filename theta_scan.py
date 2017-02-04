import sage.cboot as cb
from sage.misc.cachefunc import cached_function
from subprocess import Popen, PIPE
import re
import numpy as np
import sys
import os


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
    sn = np.cos(theta)
    cs = np.sin(theta)
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

if __name__=="__main__":
    name        = sys.argv[1]
    Lambda      = int(sys.argv[2])  # 11
    lmax        = int(sys.argv[3])  # 20
    nu_max      = int(sys.argv[4])  # 8
    precision   = int(sys.argv[5])  # 400
    resolution  = int(sys.argv[6])  # 3
    sdpbparams.append("--precision={}".format(precision))

    if len(sys.argv) > 7:
        distance = (float(sys.argv[7]), 10*float(sys.argv[7]))
    if len(sys.argv) > 8:
        ppn = int(sys.argv[8])
        sdpbparams.append("--maxThreads={}".format(ppn))
    else: distance = (0.002, 0.02)
    print "using Lambda = {}, lmax = {}, nu_max = {}, precision = {}".format(\
            Lambda, lmax, nu_max, precision)
    print "with resolution = {}, window = ({}, {}), threads = {}.".format(\
            resolution, distance[0], distance[1], ppn)

    context=cb.context_for_scalar(epsilon=0.5,Lambda=Lambda)

    Dsig = 0.518154
    Deps = 1.41267

    for delta_s in np.linspace(Dsig - distance[0], Dsig + distance[0], num = resolution):
        for delta_e in np.linspace(Deps - distance[1], Deps + distance[1], num = resolution):
            check((delta_s, delta_e), 0.969260330903202)
