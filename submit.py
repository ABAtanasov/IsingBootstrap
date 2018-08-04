# --------------------------------------------------------
# submit.py
#
# This is the module for submitting jobs to a cluster
# Written for the Yale grace HPC cluster using SLURM
# For other HPC clusters, write_jobfile must be modified
# --------------------------------------------------------

from __future__ import absolute_import

import os
import errno
import sys
import parser_tools
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
    mem       = job_params['mem']
    ndays     = job_params['ndays']
    queue     = job_params['queue']

    # Main command:
    cmd = "sage mixed_ising.py -N={} -L={} -l={} -nu={} -p={} --maxIters={} ".format(\
            name, Lambda, lmax, nu_max, precision, maxIters)\
            + "--threads={} ".format(threads)\
            + "--in_file={} --out_file=True ".format(in_file)

    # Various boolean options:
    if job_params['print_sdpb']:
        cmd += "--print_sdpb=True "
    if job_params['keepxml']:
        cmd += "--keepxml=True "
    if job_params['profile']:
        cmd += "--profile=True "
    if job_params['envelope']:
        cmd += "--envelope=True "
    if job_params['odd_scalar_gap'] is not None:
        cmd += "--odd_scalar_gap={} ".format(job_params['odd_scalar_gap'])
    if job_params['even_scalar_gap'] is not None:
        cmd += "--even_scalar_gap={} ".format(job_params['even_scalar_gap'])
    if job_params['spin_2_gap'] is not None:
        cmd += "--spin_2_gap={} ".format(job_params['spin_2_gap'])
    if job_params['max_bisections'] is not None:
        cmd += "--max_bisections={} ".format(job_params['max_bisections'])
        cmd += "--sig_spacing={} ".format(job_params['sig_spacing'])
        cmd += "--eps_spacing={} ".format(job_params['eps_spacing'])
        cmd += "--side={} ".format(job_params['side'])
    if job_params['theta_spacing'] is not None:
        cmd += "--theta_spacing={} ".format(job_params['theta_spacing'])

    mainpath = os.path.dirname(os.path.abspath(__file__))
    bsubpath = os.path.join(mainpath, "bash_scripts")
    scratchpath = os.path.join(mainpath, "scratch")
    outpath = os.path.join(mainpath, "out_files")

    mkdir_p(scratchpath)
    mkdir_p(outpath)
    mkdir_p(bsubpath)

    paths = {"main": mainpath, "sh_scripts": bsubpath,
             "scratch": scratchpath, "out": outpath}

    jobfile = write_jobfile(cmd, name, paths,
                            mem=mem, ndays=ndays, threads=threads, queue=queue)

    os.system("sbatch < {} -A poland".format(jobfile))

    print "submitted:\n  {}".format(cmd)


# --------------------------------------------------------
# This writes the entire submit .sh file
# and returns a link to that file inside bash_scripts/
#
# WARNING: Written for the slurm job scheduler
# This code would need to support other schedulers
# --------------------------------------------------------
def write_jobfile(cmd, jobname, paths,
                  nodes=10, threads=1, gpus=0, mem=8, ndays=1, queue='day'):

    days = '{}'.format(int(ndays))
    if len(days) == 1:
        days = '0' + days

    jobfile = os.path.join(paths["sh_scripts"], jobname + '.sh')

    with open(jobfile, 'w') as f:
        f.write('#! /bin/bash\n'
                + '\n'
                + '#SBATCH --job-name={}\n'.format(jobname)
                + '#SBATCH --partition={}\n'.format(queue)
                + '#SBATCH --mem-per-cpu={}\n'.format(int(mem)*1000)
                + '#SBATCH --time={}-00:00:00\n'.format(days)
                + '#SBATCH --cpus-per-task={}\n'.format(threads)
                + '#SBATCH --ntasks=1\n'
                + '#SBATCH --output={}/{}-%j.out'.format(paths['scratch'], jobname)
                + '\n'
                + 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n'
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

    # Parse lines from command line
    parser = parser_tools.build_parser()

    # Take the args as dictionary
    args = parser.parse_args().__dict__

    # Params fed into sdpb
    sdpb_params = {'Lambda': 11, 'lmax': 20, 'nu_max': 8, 'precision': 400, 'maxIters': 500, 'threads':4}

    # Params fed into the cluster/mixed_ising
    job_params = {'name':"untitled",
                  'res': [1, 1], 'theta_res':1,
                  'range': None, 'theta_range': None,
                  'dist':None, 'theta_dist':None, 'origin':None,
                  'keepxml':False, 'print_sdpb':False, 'profile':False, 'envelope':False,
                  'odd_scalar_gap': 3, 'even_scalar_gap': 3, 'spin_2_gap': None,
                  'max_bisections': None, 'sig_spacing': None, 'eps_spacing': None,
                  'theta_spacing': None, 'side': "exterior",
                  'file':None,
                  'mem':8, 'ndays':1, 'threads':4, 'queue':'pi_poland'} # This last option is cluster-dependent

    for key in sdpb_params.keys():
        if args[key] is not None:
            sdpb_params[key] = args[key]
        elif key in ["Lambda", "lmax", "nu_max"]:
            print "Warning, {} not specified. Using {} = {}.".format(
                    key, key, sdpb_params[key])

    for key in job_params.keys():
        if args[key] is not None:
            job_params[key]=args[key]

    if job_params['max_bisections']:
        assert job_params['sig_spacing'] is not None
        assert job_params['eps_spacing'] is not None

    # In the case we are only submitting one job
    if not args['batches']:
            generate_to_file(job_params)
            job_params["name"] = "{}_1of1".format(args['name'])
            job_params["file"] = "scratch/{}.pts".format(job_params['name'])
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
        job_params['file'] = "scratch/{}.pts".format(job_params['name'])
        submit_job(job_params, sdpb_params)
        print "submitted job {}".format(batch + 1)

