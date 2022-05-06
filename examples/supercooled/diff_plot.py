#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

def draw(fileName):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    jar.pop(0)

    x = []; y = []; z = []; p = []
    for line in jar:
        x.append(float(line.split(',')[0]))
        y.append(float(line.split(',')[1]))
        z.append(float(line.split(',')[2]))
        p.append(float(line.split(',')[3]))

    overT = [1/i for i in x]
    m1, b1 = np.polyfit(overT, np.log(y), 1)
    co1 = [m1*i+b1 for i in overT]
    m2, b2 = np.polyfit(overT, np.log(z), 1)
    co2 = [m2*i+b2 for i in overT]
    m3, b3 = np.polyfit(overT, np.log(p), 1)
    co3 = [m3*i+b3 for i in overT]
    co4 = [0.2*i+0.8*j for i, j in zip(co1,co2)]

    #plt.scatter(overT, y, color=plt.cm.turbo(0.3),
    #        marker='o', label=r'$D_{AA}$')
    #plt.plot(overT, np.exp(co1), color=plt.cm.turbo(0.3), ls='-')
    #plt.scatter(overT, z, color=plt.cm.turbo(0.6),
    #        marker='^', label=r'$D_{BB}$')
    #plt.plot(overT, np.exp(co2), color=plt.cm.turbo(0.6), ls='-.')
    plt.scatter(overT, p, color=plt.cm.turbo(0.9),
            marker='s', label=r'$D_{AB}$')
    plt.plot(overT, np.exp(co3), color=plt.cm.turbo(0.9), ls='-')
    plt.plot(overT, np.exp(co4), color=plt.cm.turbo(0.7), ls=':',
            label='phenom. rel.')

files = sys.argv[1:]

for file in files:
    print(file)
    draw(file)

plt.yscale('log')
plt.xlabel('1/T')
plt.ylabel(r'Diffusivity')
plt.legend()
plt.savefig('trial.pdf')
plt.show()
