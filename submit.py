from __future__ import absolute_import

import os
import errno
import sys

# Make a directory if it is not already made
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path): pass
        else: raise


def submit_job(job_params, sdpb_params):
    name   = job_params['name']
    Lambda = sdpb_params['Lambda']
    lmax   = sdpb_params['lmax']
    nu_max = sdpb_params['nu_max']
    precision = sdpb_params.get('precision', 400)
    resolution = job_params['resolution']
    distance   = job_params['distance']
    mem        = job_params.get('mem',   8)
    ndays      = job_params.get('ndays', 1)
    ppn        = job_params['ppn']

    cmd = "sage mixed_ising.py {} {} {} {} {} {} {} {}".format(\
            name, Lambda, lmax, nu_max, precision, resolution, distance, ppn)


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
            mem = mem, ndays = ndays, ppn = ppn)

    os.system("bsub < " + jobfile)
    print "submitted job {}.".format(name)

def write_jobfile(cmd, jobname, paths,
                  nodes=10, ppn=1, gpus=0, mem=8, ndays=1, queue='shared'):

    if ppn > 1:
        threads = 'export OMP_NUM_THREADS={}\n'.format(ppn)
    else:
        threads = ''

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
    if len(sys.argv) < 5:
        print "usage:"
        print "python submit.py jobname Lambda lmax nu_max precision resolution distance mem ndays ppn"
        exit(0)
    argv = sys.argv

    job_params = {'name':sys.argv[1], 'resolution':sys.argv[6], 'distance':0.002, 'mem':10, 'ndays':1, 'ppn':1}
    if len(sys.argv) > 7: job_params['distance'] = sys.argv[7]
    if len(sys.argv) > 8: job_params['mem'] = sys.argv[8]
    if len(sys.argv) > 9: job_params['ndays'] = sys.argv[9]
    if len(sys.argv) > 10: job_params['ppn'] = sys.argv[10]

    sdpb_params = {'Lambda':argv[2], 'lmax':argv[3], 'nu_max':argv[4], 'precision':argv[5]}

    print job_params
    print sdpb_params

    submit_job(job_params, sdpb_params)

