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
def print_out(string, f=None):
    if f is not None:
        f.write(string)
        f.write("\n")
    else:
        print string


def print_err(message, dumpfile, f=None):
    print_out("------------------------------------", f=f)
    print_out("An Error has occurred:", f=f)
    print_out(message, f=f)
    print_out("------------------------------------", f=f)
    print_out(dumpfile, f=f)


def make_decimal(number):
    return "%.3E" % Decimal(number)


def print_point(deltas, out, name, durations, profile=True, f=None):
    if len(deltas) == 3:
        message = "({}, {}, {}) ".format(deltas[0], deltas[1], deltas[2])
    else:
        message = "({}, {}) ".format(deltas[0], deltas[1])

    sol = re.compile(r'found ([^ ]+) feasible').search(out)
    if not sol:
        print_err("The out file wasn't right: ", out, f=f)
        os.system("rm scratch/{}.ck".format(name))
        message += "is not excluded. "
        print_out(message, f=f)
        write_update(f)
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
            print_err("No MPI speedup message", out, f=f)
        else:
            speedup = speedup.groups()[0]
            message += "| MPI speedup = {}% ".format(speedup)
        if not dual or not primal:
            print_err("No line for either dual or primal error", out, f=f)
        else:
            dual = dual.groups()[0]
            primal = primal.groups()[0]
            message += "| Errors: Dual = {}, Primal = {} ".format(
                    make_decimal(dual), make_decimal(primal))

    print_out(message, f=f)
    write_update(f)

    return excluded
