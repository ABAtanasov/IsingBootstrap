# Imports a job

from __future__ import absolute_import

import os
import errno
import sys
import argparse
import numpy as np


# --------------------------------------------------------
# Makes a numpy array from a to b with a given resolution
# --------------------------------------------------------
def mkrange(a, b, resolution):
    if resolution == 1:
        return np.array([0.5*(a+b)])
    else:
        return np.linspace(a, b, num = resolution)

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

    range     = job_params['range']
    theta_range = job_params['theta_range']
    res       = job_params['res']
    theta_res = job_params['theta_res']
    keepxml   = job_params['keepxml']
    printxml  = job_params['printxml']
    mem       = job_params['mem']
    ndays     = job_params['ndays']
    queue     = job_params['queue']

    scaling_info = theta_info = " "
    if range:
        scaling_info = "--range {} {} {} {} --res {} {} ".format(\
                range[0], range[1], range[2], range[3], res[0], res[1])
    if theta_range:
        theta_info = "--theta_range {} {} --theta_res={} ".format(\
                theta_range[0], theta_range[1], theta_res)

    cmd = "sage mixed_ising.py -N={} -L={} -l={} -nu={} -p={} --keepxml={} --printxml={} --maxIters={} ".format(\
            name, Lambda, lmax, nu_max, precision, keepxml, printxml, maxIters)\
            + scaling_info\
            + theta_info\
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
            + '{} > {}/{}.out 2>&1\n'.format(cmd, paths["out"], jobname)
            + '\n'
            + 'exit 0;\n'
            )

    return jobfile

    sdpb_params = {'Lambda':11, 'lmax':20, 'nu_max':8, 'precision':400, 'maxIters':500, 'threads':4}

    # Params fed into the cluster
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
    parser.add_argument("-N", "--name", type = str,\
            help="name for the associated files")
    parser.add_argument("-L","--Lambda", type = int, \
           help="maximum derivative order")
    parser.add_argument("-l", "--lmax", type = int, \
            help="angular momentum cutoff")
    parser.add_argument("-nu", "--nu_max", type = int, \
            help="maximum number of poles")
    parser.add_argument("--res", type = int, nargs = 2,\
            help="number of sampling points along each axis")
    parser.add_argument("--theta_res", type = int,\
            help="number of sampling points over use_theta")
    parser.add_argument("--dist", type = float,\
            help="distance of Delta_sigma window from the 3D Ising point")
    parser.add_argument("--theta_dist", type = float,\
            help="distance of use_theta window from the 3D Ising use_theta")
    parser.add_argument("--range", type = float, nargs = 4,\
            help="4 floats xmin xmax ymin ymax")
    parser.add_argument("--origin", type = float, nargs = 2,\
            help="2 floats x_origin y_origin")
    parser.add_argument("--theta_range", type = float, nargs = 2,\
            help="2 floats theta_min theta_max")
    parser.add_argument("--keepxml", type=bool, \
                        help="Do we keep the xml? Default is no.")
    parser.add_argument("--printxml", type=bool, \
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

    Dsig = 0.518154
    Deps = 1.41267
    distance   = (0.002, 0.02)
    theta0 = 0.969260330903202

    if args['origin'] and not args['range']:
        Dsig, Deps = args['origin']
    if args['dist'] and not args['range']:
        dist = float(args['dist'])
        distance = (dist, 10*dist)
    sig_min = Dsig - distance[0]
    sig_max = Dsig + distance[0]
    eps_min = Deps - distance[1]
    eps_max = Deps + distance[1]
    if args['range']:
        sig_min, sig_max, eps_min, eps_max = args['range']
    scaling_range = [sig_min, sig_max, eps_min, eps_max]

    use_theta = False
    if args['theta_range']:
        use_theta = True
        theta_min, theta_max = args['theta_range']
    elif args['theta_dist']:
        use_theta = True
        theta_min = theta0 - args['theta_dist']
        theta_max = theta0 + args['theta_dist']
    theta_range = None
    if use_theta:
        theta_range = [theta_min, theta_max]


    # Params fed into sdpb
    sdpb_params = {'Lambda': 11, 'lmax': 20, 'nu_max': 8, 'precision': 400, 'maxIters': 500, 'threads':4}

    # Params fed into the cluster/mixed_ising
    job_params = {'name':"untitled",
            'res':[1, 1], 'theta_res':1,
            'range': None, 'theta_range': None,
            'keepxml':False, 'printxml':False,
            'mem':8, 'ndays':1,'queue':'shared'} # This last option is cluster-dependent

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
            print "Using", sdpb_params
            print "with", job_params
            submit_job(job_params, sdpb_params)
            exit(0)

    # In the case where we are submitting multiple jobs
    batches = args['batches']
    res = job_params['res']
    theta_res = job_params['theta_res']

    print "Using", sdpb_params

    if len(batches) == 1:
        theta_bits = mkrange(theta_range[0], theta_range[1], batches+1)
        theta_ranges = [[theta_bits[i+1], theta_bits[i]] for i in range(batches)]
        for batch in range(batches):
            job_params['name'] = "{}_{}-{}of{}".format(
                                args["name"], batch + 1, batches)
            job_params['theta_range'] = [theta_ranges[batch][0], theta_ranges[batch][1]]
            job_params['theta_res']   = [theta_res/batches]
            additional_info = " --res {} {}".format(res[0], res[1]) \
                    +" --theta_res {}".format(theta_res/batches)\
                    +" --range {} {} {} {}".format(\
                    scaling_range[0], scaling_range[1], scaling_range[2], scaling_range[3])\
                    +" --theta_range {} {}".format(\
                    theta_ranges[batch][0], theta_ranges[batch][1])\
                    +" -N {}_{}of{}".format(args["name"], batches)
            os.system(cmd + additional_info)
            print "submitted:\n {}".format(cmd + additional_info)

    elif len(batches) == 2:
        total_batches = batches[0] * batches[1]
        sig_bits = mkrange(scaling_range[0], scaling_range[1], batches[0]+1)
        eps_bits = mkrange(scaling_range[2], scaling_range[3], batches[1]+1)
        sig_ranges = [[sig_bits[i], sig_bits[i+1]] for i in range(batches[0])]
        eps_ranges = [[eps_bits[i], eps_bits[i+1]] for i in range(batches[1])]
        for sig_batch in range(batches[0]):
            for eps_batch in range(batches[1]):
                job_params['name'] = "{}_{}-{}of{}".format(
                                args["name"], sig_batch+1, eps_batch+1, total_batches)
                job_params['range'] = [sig_ranges[sig_batch][0], sig_ranges[sig_batch][1],
                                       eps_ranges[eps_batch][0], eps_ranges[eps_batch][1]]
                job_params['res']   = [res[0]/batches[0], res[1]/batches[1]]
                print "submitted job with", job_params
                submit_job(job_params, sdpb_params)