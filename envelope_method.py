# -----------------------------------------------------------------
# envelope_method.py
#
# Implementing the envelope method discussed previously with D. Poland and A. Hillman
# to significantly decrease the number of points that must be checked in future theta scans
#
# -----------------------------------------------------------------

from point_generator import array2dict
import random

def approx(a, b):
    return abs(a - b) < 1e-10

def get_step(theta_fiber):
    assert sorted(theta_fiber)
    if len(theta_fiber) == 1:
        return None
    else:
        step = min([t2 - t1 for t1, t2 in zip(theta_fiber[:-1], theta_fiber[1:])])
        if step == 0.0:
            raise RuntimeError("There are repeated thetas in the list.")
        else:
            return step


def get_chunks(theta_fiber):
    step = get_step(theta_fiber)
    chunks = []
    chunk = [theta_fiber[0]]
    if step is None:
        return [chunk]
    for i in range(1, len(theta_fiber)):
        if not approx(theta_fiber[i] - theta_fiber[i-1], step):
            chunks.append(chunk)
            chunk = [theta_fiber[i]]
        else:
            chunk.append(theta_fiber[i])
    chunks.append(chunk)
    return chunks


def erode(check, base_point, chunk, f=None):
    assert sorted(chunk)

    # Go down from the top until something isn't excluded
    for above in range(len(chunk)-1, -1, -1):
        if not check(base_point, theta=chunk[above], f=f):
            break

    # Go down until something isn't excluded
    for below in range(0, above):
        if not check(base_point, theta=chunk[below], f=f):
            break


def envelope_loop(check, points, f=None):
    base_points = array2dict(points)
    for base_point, thetas in base_points.iteritems():
        # Over each base point, you can assume positive values of theta (for now.. with this batch..)
        # with uniform step size (in fact 0.01 for this batch)
        # From there, get the connected components and erode
        thetas.sort()
        chunks = get_chunks(thetas)
        for chunk in chunks:
            erode(check, base_point, chunk, f=f)
