#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

files = sys.argv[1:]

for file in files:
    temp = float(file[3:6])
    with open(file, 'r') as f:
        jar = f.readlines()
    jar = jar[1:-1]

    x = []; y = []
    for line in jar:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

    tint = plt.cm.turbo((temp-0.7)/0.7*0.8+0.1)
    plt.plot(x, y, label=f'{temp}', color=tint)

plt.xlabel('Time')
plt.ylabel('Diffusion coefficient')
plt.legend()
plt.savefig('diffCoeff.pdf')
plt.show()
