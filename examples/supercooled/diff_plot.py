#!/usr/bin/env python3

import sys
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
    #plt.plot(overT, y, color=plt.cm.turbo(0.3),
    #        marker='o', label=r'$D_{AA}$', ls='--')
    #plt.plot(overT, z, color=plt.cm.turbo(0.6),
    #        marker='^', label=r'$D_{BB}$', ls='-.')
    plt.plot(overT, p, color=plt.cm.turbo(0.9),
            marker='s', label=r'$D_{AB}$', ls='-')

files = sys.argv[1:]

for file in files:
    print(file)
    draw(file)

plt.xlabel('1/T')
plt.ylabel(r'Diffusivity')
plt.legend()
plt.savefig('trial.pdf')
plt.show()
