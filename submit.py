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
# Given the job params and sdpb params this will set up
# everything for the cluster to run cboot and sdpb
# --------------------------------------------------------
def submit_job(job_params, sdpb_params):
    name      = job_params['name']
    program   = job_params['program']
    Lambda    = sdpb_params['Lambda']
    lmax      = sdpb_params['lmax']
    nu_max    = sdpb_params['nu_max']
    precision = sdpb_params['precision']
    res       = job_params['res']
    dist      = job_params['dist']
    theta_res = job_params['theta_res']
    theta_dist= job_params['theta_dist']
    range     = job_params['range']
    theta_range = job_params['theta_range']
    mem       = job_params['mem']
    ndays     = job_params['ndays']
    threads   = job_params['threads']
    queue     = job_params['queue']

    if range:
        scaling_info = "--range {} {} {} {} --res {} {} ".format(\
                range[0], range[1], range[2], range[3], res[0], res[1])
    else:
        scaling_info = "--dist={} --res {} {} ".format(dist, res[0], res[1])

    if theta_range:
        theta_info = "--theta_range {} {} --theta_res={} ".format(\
                theta_range[0], theta_range[1], theta_res)
    else:
        theta_info = "--theta_dist {} --theta_res={} ".format(\
                theta_dist, theta_res)

    if program == "mixed_ising":
        theta_info = ""

    cmd = "sage {}.py -N={} -L={} -l={} -nu={} -p={} ".format(\
            program, name, Lambda, lmax, nu_max, precision)\
            + scaling_info\
            + theta_info\
            + "--threads={}".format(threads)


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
            + '{} > {}/{}.out 2>&1\n'.format(cmd, paths["out"], jobname)
            + '\n'
            + 'exit 0;\n'
            )

    return jobfile

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

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
    parser.add_argument("--theta_res", type = int,\
            help="number of sampling points over theta")
    parser.add_argument("--dist", type = float,\
            help="distance of Delta_sigma window from the 3D Ising point")
    parser.add_argument("--theta_dist", type = float,\
            help="distance of theta window from the 3D Ising theta")
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

    # If no flags are given, print the help menu instead:
    if len(sys.argv) == 1:
        os.system("python {} -h".format(os.path.abspath(__file__)))
        exit(0)

    distance   = (0.002, 0.02)
    if args['range']:
        sig_min = args['range'][0]
        sig_max = args['range'][1]
        eps_min = args['range'][2]
        eps_max = args['range'][3]
    elif args['dist']:
        dist = float(args['dist'])
        distance = (dist, 10*dist)

    # Params fed into sdpb
    sdpb_params = {'Lambda':11, 'lmax':20, 'nu_max':8, 'precision':400}

    # Params fed into the cluster
    job_params = {'name':"untitled", 'program':"mixed_ising",\
            'res':[1, 1], 'theta_res':1, 'dist':0.002,\
            'theta_dist':0.1, 'mem':8, 'ndays':1, 'threads':4,\
            'range':None, 'theta_range':None,\
            'queue':'shared'}

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

