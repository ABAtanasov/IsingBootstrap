# --------------------------------------------------------
# printing_tools.py
#
# This is the module for printing various messages
# to various output files f_out
# (or STDOUT if f_out is None)
# --------------------------------------------------------
import time
import re
import os
from decimal import Decimal

# --------------------------------------------------------
# Updates a file being written to
# --------------------------------------------------------
def write_update(f_out):
    if f_out is not None:
        f_out.flush()
        time.sleep(4)

# --------------------------------------------------------
# Prints to the file f or stdout if no f is specified
# --------------------------------------------------------
def print_out(string, f_out=None):
    if f_out is not None:
        f_out.write(string)
        f_out.write("\n")
    else:
        print string


def print_err(message, dumpfile, f_out=None):
    print_out("------------------------------------", f_out=f_out)
    print_out("An Error has occurred:", f_out=f_out)
    print_out(message, f_out=f_out)
    print_out("------------------------------------", f_out=f_out)
    print_out(dumpfile, f_out=f_out)


def make_decimal(number):
    return "%.3E" % Decimal(number)


def print_point(deltas, out, name, durations, profile=True, f_out=None):
    if len(deltas) == 3:
        message = "({}, {}, {}) ".format(deltas[0], deltas[1], deltas[2])
    else:
        message = "({}, {}) ".format(deltas[0], deltas[1])

    sol = re.compile(r'found ([^ ]+) feasible').search(out)
    if not sol:
        print_err("The out file wasn't right: ", out, f_out=f_out)
        os.system("rm scratch/{}.ck".format(name))
        message += "is not excluded. "
        print_out(message, f_out=f_out)
        write_update(f_out)
        return False
    sol = sol.groups()[0]
    if sol == "dual":
        message += "is excluded. "
        excluded = True
    elif sol == "primal":
        message += "is not excluded. "
        excluded = False
    else:
        raise RuntimeError

    if profile:
        speedup = re.compile(r'Solver runtime.+ CPU \(([^ ]+)%\)').search(out)
        dual = re.compile(r'dualError[ ]+= ([^\s]+)').search(out)
        primal = re.compile(r'primalError[ ]+= ([^\s]+)').search(out)
        message += "| cboot = {}s, sdpb = {}s ".format(durations[0], durations[1])
        if not speedup:
            print_err("No MPI speedup message", out, f_out=f_out)
        else:
            speedup = speedup.groups()[0]
            message += "| MPI speedup = {}% ".format(speedup)
        if not dual or not primal:
            print_err("No line for either dual or primal error", out, f_out=f_out)
        else:
            dual = dual.groups()[0]
            primal = primal.groups()[0]
            message += "| Errors: Dual = {}, Primal = {} ".format(
                    make_decimal(dual), make_decimal(primal))

    print_out(message, f_out=f_out)
    write_update(f_out)

    return excluded
