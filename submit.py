# Imports a job

from __future__ import absolute_import

import os
import errno
import sys
import argparse
import numpy as np
from point_generator import generate_to_file


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
# Given the job params and sdpb params this will set up
# everything for the cluster to run cboot and sdpb
# --------------------------------------------------------
def submit_job(job_params, sdpb_params):
    name      = job_params['name']

    Lambda    = sdpb_params['Lambda']
    lmax      = sdpb_params['lmax']
    nu_max    = sdpb_params['nu_max']
    precision = sdpb_params['precision']
    maxIters  = sdpb_params['maxIters']
    threads   = job_params['threads']
    in_file   = job_params['file']
    keepxml   = job_params['keepxml']
    print_sdpb  = job_params['print_sdpb']
    mem       = job_params['mem']
    ndays     = job_params['ndays']
    queue     = job_params['queue']

    cmd = "sage mixed_ising.py -N={} -L={} -l={} -nu={} -p={} --keepxml={} --print_sdpb={} --maxIters={} --in_file={} ".format(\
            name, Lambda, lmax, nu_max, precision, keepxml, print_sdpb, maxIters, in_file)\
            + "--threads={} ".format(threads)


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
            mem = mem, ndays = ndays, threads = threads, queue = queue)

    os.system("bsub < " + jobfile)

    print "submitted:\n  {}".format(cmd)


# --------------------------------------------------------
# This writes the entire submit .sh file
# and returns a link to that file inside bash_scripts/
# --------------------------------------------------------
def write_jobfile(cmd, jobname, paths,
                  nodes=10, threads=1, gpus=0, mem=8, ndays=1, queue='shared'):

    exp_threads = 'export OMP_NUM_THREADS={}\n'.format(threads)

    jobfile = os.path.join(paths["sh_scripts"], jobname + '.sh')

    with open(jobfile, 'w') as f:
        f.write('#! /bin/bash\n'
            + '\n'
            #+ '#BSUB -l nodes={}:ppn={}\n'.format(nodes, ppn)
            + '#BSUB -M {}\n'.format(int(mem)*1000)          #good
            + '#BSUB -W {}:00\n'.format(24*int(ndays))       #good
            + '#BSUB -n {}\n'.format(threads)
            + '#BSUB -q {}\n'.format(queue)
            + '#BSUB -J {}\n'.format(jobname)                #good
            + '#BSUB -e {}/{}.e%J\n'.format(paths["scratch"], jobname)
            + '#BSUB -o {}/{}.o%J\n'.format(paths["scratch"], jobname)
            + '\n'
            + exp_threads                                    #good
            + '. ~/.bashrc\n'
            + '{} 2>&1\n'.format(cmd)
            + '\n'
            + 'exit 0;\n'
            )

    return jobfile


if __name__ == "__main__":

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("python {} -h".format(os.path.abspath(__file__)))
        exit(0)

    parser = argparse.ArgumentParser()

    # --------------------------------------
    # Args for submit.py
    # --------------------------------------
    parser.add_argument("-B", "--batches", type=int, nargs='+', \
            help="info for how jobs are submitted into batches")

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
    parser.add_argument("-f", "--file", type=str,
                        help="optional filename from which we read the points")
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
    parser.add_argument("--keepxml", type=bool,
                        help="Do we keep the xml? Default is no.")
    parser.add_argument("--print_sdpb", type=bool,
                        help="Do we print the sdpb output? Default is no.")

    # --------------------------------------
    # Args for sdpb
    # --------------------------------------
    parser.add_argument("-p", "--precision", type =int, \
            help="working precision for sdpb calculations")
    parser.add_argument("-i", "--maxIters", type=int, \
            help="max number of sdpb iterations")
    parser.add_argument("--threads", type=int, \
                        help="maximum threads used by OpenMP")

    # --------------------------------------
    # Args for the cluster
    # --------------------------------------
    parser.add_argument("--mem", type = int,\
            help="maximum memory allocated per node in cluster")
    parser.add_argument("--ndays", type = int,\
            help="number of days to run process on cluster")
    parser.add_argument("-q", "--queue", type = str,\
            help="queue to submit to")

    # Take the args as dictionary
    args = parser.parse_args().__dict__

    # Params fed into sdpb
    sdpb_params = {'Lambda': 11, 'lmax': 20, 'nu_max': 8, 'precision': 400, 'maxIters': 500, 'threads':4}

    # Params fed into the cluster/mixed_ising
    job_params = {'name':"untitled",
            'res':[1, 1], 'theta_res':1,
            'range': None, 'theta_range': None,
            'dist':None, 'theta_dist':None, 'origin':None,
            'keepxml':False, 'print_sdpb':False, 'file':None,
            'mem':8, 'ndays':1,'queue':'shared', 'threads':4} # This last option is cluster-dependent

    for key in sdpb_params.keys():
        if args[key]:
            sdpb_params[key] = args[key]
        elif key in ["Lambda", "lmax", "nu_max"]:
            print "Warning, {} not specified. Using {} = {}.".format(
                    key, key, sdpb_params[key])

    for key in job_params.keys():
        if args[key]:
            job_params[key]=args[key]

    # In the case we are only submitting one job
    if not args['batches']:
            generate_to_file(job_params)
            submit_job(job_params, sdpb_params)
            exit(0)

    # In the case where we are submitting multiple jobs

    batches = args['batches']
    f_in = None
    if job_params["file"] is not None:
        f_in = open(job_params["file"], "r")

    generate_to_file(job_params, batches=batches, f_in=f_in)

    for batch in range(batches):
        job_params['name'] = "{}_{}of{}".format(args["name"], batch + 1, batches)
        job_params['file'] = "scratch/{}".format(job_params['name'])
        submit_job(job_params, sdpb_params)
        print "submitted job {}".format(batch + 1)
