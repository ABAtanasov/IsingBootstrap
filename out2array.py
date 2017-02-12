import sys
import re

lines = sys.stdin.readlines()

inclusion = re.compile("is not excluded")
number_data = re.compile("[\d]+.[\d]+")
thetas = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'thetas':
        thetas = True


sigmas = []
epsilons = []
thetas = []

for line in lines:
    if inclusion.search(line) is not None:
        data = number_data.findall(line)
        sigmas.append(float(data[0]))
        epsilons.append(float(data[1]))
        if thetas: thetas.append(float(data[2]))


print sigmas
print epsilons
if thetas: print thetas
