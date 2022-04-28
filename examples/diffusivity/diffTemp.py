#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

files = sys.argv[1:]

temps = []
diffs = []
for file in files:
    temp = float(file[3:6])
    temps.append(temp)
    with open(file, 'r') as f:
        jar = f.readlines()
    diffs.append(float(jar[-1].split()[1]))

overT = [1/T for T in temps]
k = 5.732e-8
diffuse = [k*x for x in diffs]
plt.scatter(overT, diffuse, color=plt.cm.turbo(0.2))

m, b = np.polyfit(overT, np.log(diffuse), 1)
ordinates = [m*x+b for x in overT]
plt.plot(overT, np.exp(ordinates), color=plt.cm.turbo(0.3))

plt.yscale('log')
plt.xlabel('1/T')
plt.ylabel('Diffusion coefficient')
plt.tight_layout()
plt.savefig('overT.pdf')
