

# Performs a round of bisection for a given point
def bisect(check, point, spacing, num_bisections=0, max_bisections=5, f=None, side="exterior"):
    if num_bisections >= max_bisections:
        return

    excluded = check(point, f=f)
    keep_going = True
    if side == "exterior" and not excluded:
        keep_going = False
    if side == "interior" and excluded:
        keep_going = False

    if keep_going and num_bisections == 0:
        new_spacing = spacing
        new_num_bisections = num_bisections
    else:
        new_spacing = [diff/2.0 for diff in spacing]
        new_num_bisections = num_bisections + 1

    if keep_going:
        new_point = [pt + diff for pt, diff in zip(point, new_spacing)]
    else:
        new_point = [pt - diff for pt, diff in zip(point, new_spacing)]

    bisect(check, new_point, new_spacing,
           num_bisections=new_num_bisections, max_bisections=max_bisections,
           f=f, side=side)


# Loops over an input set of points, performing the bisection method individually on each one
def bisect_loop(check, points, spacing, max_bisections, f=None, side="exterior"):
    assert side == "exterior" or side == "interior"
    for point in points:
        bisect(check, point, spacing=spacing, max_bisections=max_bisections, f=f, side=side)

