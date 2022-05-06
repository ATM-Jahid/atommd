#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

masses = []
A = []; B = []; AB = []
M = []

def extr(fileName):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    jar.pop(0)

    for line in jar:
        tmp = line.split(',')
        mass = float(tmp[0])
        masses.append(mass)
        Q = 1/0.2/0.8 * (mass*0.8+0.2)**2
        A.append(float(tmp[1]))
        B.append(float(tmp[2]))
        AB.append(float(tmp[3]))
        M.append(0.2*float(tmp[1])+0.8*float(tmp[2]))

files = sys.argv[1:]
for file in files:
    extr(file)

R = [x/y for x, y in zip(A, B)]

for a, b, c in zip(masses, R, AB):
    print(f'{a} & {b} & {c}')
plt.plot(masses, R, marker='s', ls='-.', label=r'$D_A/D_B$')
plt.plot(masses, AB, marker='P', label=r'$D_{AB}$')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mass ratio')
plt.ylabel('Diffusion coefficients')
plt.legend()
plt.savefig('trial.pdf')
plt.show()
