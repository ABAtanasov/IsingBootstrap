import sys
import re

lines = sys.stdin.readlines()

inclusion = re.compile("is not excluded")
number_data = re.compile("[\d]+.[\d]+")


sigmas = []
epsilons = []
thetas = []

for line in lines:
    if inclusion.search(line) is not None:
        data = number_data.findall(line)
        sigmas.append(float(data[0]))
        epsilons.append(float(data[1]))
        thetas.append(float(data[2]))


print sigmas
print epsilons
print thetas
