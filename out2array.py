import sys
import re

lines = sys.stdin.readlines()

inclusion = re.compile("is not excluded")
number_data = re.compile("[\d]+.[\d]+")
theta = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'theta':
        theta = True


sigmas = []
epsilons = []
thetas = []

for line in lines:
    if inclusion.search(line) is not None:
        data = number_data.findall(line)
        if len(data) > 2: theta = True
        sigmas.append(float(data[0]))
        epsilons.append(float(data[1]))
        if theta: thetas.append(float(data[2]))


print sigmas
print epsilons
if theta: print thetas
