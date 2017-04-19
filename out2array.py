import sys
import re

lines = sys.stdin.readlines()

inclusion = re.compile("is not excluded")
number_data = re.compile("[\d]+.[\d]+")
use_theta = False

sigmas = []
epsilons = []
thetas = []

for line in lines:
    if inclusion.search(line) is not None:
        data = number_data.findall(line)
        if len(data) > 2: use_theta = True
        sigmas.append(float(data[0]))
        epsilons.append(float(data[1]))
        if use_theta:
            thetas.append(float(data[2]))


print sigmas
print epsilons
if use_theta:
    print thetas
