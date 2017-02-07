from __future__ import absolute_import

import os
import errno
import sys
import argparse


# --------------------------------------------------------
# Makes a directory if it is not already made
# --------------------------------------------------------
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path): pass
        else: raise

# --------------------------------------------------------
# Given the job params and sdpb params this will set up everything
# for the cluster to run cboot and sdpb
# --------------------------------------------------------
def submit_job(job_params, sdpb_params):
    name      = job_params['name']
    Lambda    = sdpb_params['Lambda']
    lmax      = sdpb_params['lmax']
    nu_max    = sdpb_params['nu_max']
    precision = sdpb_params['precision']
    res       = job_params['res']
    theta_res = job_params['theta_res']
    dist      = job_params['dist']
    mem       = job_params['mem']
    ndays     = job_params['ndays']
    threads   = job_params['threads']

    cmd = "sage theta_scan.py -n={} -L={} -l={} -nu={} -p={}".format(\
            name, Lambda, lmax, nu_max, precision)\
            + " --res={} --theta_res={} --dist={} --threads={}".format(\
            res, theta_res, dist, threads)


    mainpath = os.path.dirname(os.path.abspath(__file__))
    bsubpath = os.path.join(mainpath, "bash_scripts")
    scratchpath = os.path.join(mainpath, "scratch")
    outpath = os.path.join(mainpath, "out_files")

    mkdir_p(scratchpath)
    mkdir_p(outpath)
    mkdir_p(bsubpath)

    paths = {"main":mainpath, "sh_scripts":bsubpath, \
             "scratch":scratchpath, "out":outpath}

    jobfile = write_jobfile(cmd, name, paths,\
            mem = mem, ndays = ndays, threads = threads)

    os.system("bsub < " + jobfile)

    print "submitted:\n  {}.".format(cmd)


# --------------------------------------------------------
# This writes the entire submit .sh file
# and returns a link to that file inside bash_scripts/
# --------------------------------------------------------
def write_jobfile(cmd, jobname, paths,
                  nodes=10, threads=1, gpus=0, mem=8, ndays=1, queue='shared'):

    threads = 'export OMP_NUM_THREADS={}\n'.format(threads)

    jobfile = os.path.join(paths["sh_scripts"], jobname + '.sh')

    with open(jobfile, 'w') as f:
        f.write('#! /bin/bash\n'
            + '\n'
            #+ '#BSUB -l nodes={}:ppn={}\n'.format(nodes, ppn)
            + '#BSUB -M {}\n'.format(int(mem)*1000)          #good
            + '#BSUB -W {}:00\n'.format(24*int(ndays))       #good
            + '#BSUB -q {}\n'.format(queue)
            + '#BSUB -J {}\n'.format(jobname)           #good
            + '#BSUB -e {}/{}.e%J\n'.format(paths["scratch"], jobname)
            + '#BSUB -o {}/{}.o%J\n'.format(paths["scratch"], jobname)
            + '\n'
            + threads                                   #good
            + '. ~/.bashrc\n'
            + '{} > {}/{}.out 2>&1\n'.format(cmd, paths["out"], jobname)
            + '\n'
            + 'exit 0;\n'
            )

    return jobfile

if __name__ == "__main__":

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
    parser.add_argument("--theta_dist", type = float,\
            help="distance of theta window from the 3D Ising theta")
    parser.add_argument("--mem", type = int,\
            help="maximum memory allocated per node in cluster")
    parser.add_argument("--ndays", type = int,\
            help="number of days to run process on cluster")
    parser.add_argument("--threads", type = int, \
            help="maximum threads used by OpenMP")
    args = parser.parse_args().__dict__

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("python {} -h".format(os.path.abspath(__file__)))
        exit(0)

    # Params fed into sdpb
    sdpb_params = {'Lambda':11, 'lmax':20, 'nu_max':8, 'precision':400}

    # Params fed into the cluster
    job_params = {'name':"untitled", 'res':1, 'theta_res':1, 'dist':0.002,\
            'theta_dist':0.1, 'mem':8, 'ndays':1, 'threads':4}

    for key in sdpb_params.keys():
        if args[key]:
            sdpb_params[key] = args[key]
        elif key != "precision":
            print "Warning, {} not specified. Using {} = {}.".format(\
                    key, key, sdpb_params[key])

    for key in job_params.keys():
        if args[key]:
            job_params[key]=args[key]

    print "Using", sdpb_params
    print "with", job_params

    submit_job(job_params, sdpb_params)

