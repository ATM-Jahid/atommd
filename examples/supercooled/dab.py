#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

def dot_draw(fileName, itr):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    ind = fileName.find('_') + 1
    temp = float(fileName[ind:ind+3])
    jar = jar[2:]

    x = []; a = []; b = []; c= []
    for line in jar:
        x.append(float(line.split()[0]))
        a.append(float(line.split()[8]))
        b.append(float(line.split()[9]))
        c.append(float(line.split()[10]))

    d = [0.2*i + 0.8*j for i, j in zip(a, b)]
    plt.plot(x, c, color=plt.cm.turbo(itr/5),
            ls='-', label=f'{temp} '+r'$D_{AB}$')
    plt.plot(x, d, color=plt.cm.turbo(itr/5),
            ls=':', label=f'{temp} rel.')
    print(f'{a[-1]} & {b[-1]} & {d[-1]} & {c[-1]}')

files = sys.argv[1:]

for file, itr in zip(files, range(5)):
    print(file, itr)
    dot_draw(file, itr)

plt.xlabel('time')
plt.ylabel('Binary diffusion coefficient')
plt.legend()
plt.savefig('trial.pdf')
plt.show()
