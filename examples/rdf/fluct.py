#!/usr/bin/env python3

import sys
import numpy as np
from statistics import mean, stdev
import matplotlib.pyplot as plt

files = sys.argv[1:]

N = 864
temps = []
specs = []
for file in files:
    temp = float(file[3:6])
    temps.append(temp)

    with open(file, 'r') as f:
        jar = f.readlines()
    jar = jar[200:-1]

    kin = []
    for line in jar:
        kin.append(float(line.split()[3]))

    foo = stdev(kin)**2
    specs.append(3/2 * (1 - 2*N*foo/3/temp**2)**(-1))

plt.scatter(temps, specs, s=15, color=plt.cm.YlOrRd(0.9))
plt.plot(temps, specs, color=plt.cm.YlOrRd(0.4))

plt.xlabel('Temperature')
plt.ylabel('Specific heat')
plt.savefig('flucSpec.pdf')
plt.show()
