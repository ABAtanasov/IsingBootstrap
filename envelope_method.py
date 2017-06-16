# -----------------------------------------------------------------
# envelope_method.py
#
# Implementing the envelope method discussed previously with D. Poland and A. Hillman
# to significantly decrease the number of points that must be checked in future theta scans
#
# -----------------------------------------------------------------

from point_generator import array2dict
from mixed_ising import check


def envelope(points, f=None):
    """

    :param points:
    :param f:
    :return:
    """
    base_points = array2dict(points)
    for base_point, thetas in base_points.iteritems():
        # Go down from the top until something isn't excluded
        above = len(thetas)
        for above in range(len(thetas), 0, -1):
            if not check(base_point, theta=thetas[above], f=f):
                break
        # Go down until something isn't excluded
        for below in range(0, above):
            if not check(base_point, theta=thetas[below], f=f):
                break
