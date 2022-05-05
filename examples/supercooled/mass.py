#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

masses = []
A = []; B = []; AB = []

def extr(fileName):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    ind = fileName.find('_') + 1
    ind2 = fileName.find('dfs') - 1
    mass = float(fileName[ind:ind2])
    masses.append(mass)

    jar = jar[-1]
    tmp = jar.split()

    A.append(float(tmp[1]))
    B.append(float(tmp[2]))
    AB.append(float(tmp[3]))

files = sys.argv[1:]
for file in files:
    extr(file)

R = [x/y for x, y in zip(A, B)]

plt.plot(masses, R, marker='o')
plt.plot(masses, AB, marker='o')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
