# -----------------------------------------------------------------
# envelope_method.py
#
# Implementing the envelope method
# This takes an input dict of points with a base (1 or 2)-tuple and
# an array of (2 or 3)-tuples with that base tuple as their first component
#
# For each base point, we `erode away` this array of possible points
# both from the highest value down and from the lowest value up
# by excluding points using sdpb until somethnig is not excluded
# -----------------------------------------------------------------

from point_generator import array2dict2D, array2dict3D


def approx(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <  max(rel_tol * max(abs(a), abs(b)), abs_tol)


def get_step(fiber):
    assert sorted(fiber)
    if len(fiber) == 1:
        return None
    else:
        step = min([t2 - t1 for t1, t2 in zip(fiber[:-1], fiber[1:])])
        if step == 0.0:
            raise RuntimeError("There are repeated values in the list.")
        else:
            return step


def get_chunks(fiber):
    step = get_step(fiber)
    chunks = []
    chunk = [fiber[0]]
    if step is None:
        return [chunk]
    for i in range(1, len(fiber)):
        if not approx(fiber[i] - fiber[i-1], step):
            chunks.append(chunk)
            chunk = [fiber[i]]
        else:
            chunk.append(fiber[i])
    chunks.append(chunk)
    return chunks


def erode(check, base_point, chunk, f=None):
    assert sorted(chunk)

    # Go down from the top until something isn't excluded
    for above in range(len(chunk)-1, -1, -1):
        point = tuple(list(base_point) + [chunk[above]])
        if not check(point, f=f):
            break

    # Go down until something isn't excluded
    for below in range(0, above):
        point = tuple(list(base_point) + [chunk[below]])
        if not check(point, f=f):
            break


def envelope_loop3D(check, points, f=None):
    base_points = array2dict3D(points)
    for base_point, thetas in base_points.iteritems():
        # Over each base point, you can assume positive values of theta (for now.. with this batch..)
        # with uniform step size (in fact 0.01 for this batch)
        # From there, get the connected components and erode
        thetas.sort()
        chunks = get_chunks(thetas)
        for chunk in chunks:
            erode(check, base_point, chunk, f=f)


def envelope_loop2D(check, points, f=None):
    base_points = array2dict2D(points)
    for base_point, epsilons in base_points.iteritems():
        epsilons.sort()
        chunks = get_chunks(epsilons)
        for chunk in chunks:
            erode(check, base_point, chunk, f=f)

