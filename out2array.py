import sys
import re

if __name__ == "__main__":

    inclusion = re.compile("is not excluded")
    number_data = re.compile("[\d]+.[\d]+")

    points = []
    for line in sys.stdin.readlines():
        if inclusion.search(line) is not None:
            data = number_data.findall(line)
            points.append(map(lambda x: float(x), data))

    output = "{}".format(points)
    output = output.replace("[", "{")
    output = output.replace("]", "}")

    print output
