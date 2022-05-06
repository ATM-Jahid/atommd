#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

files = sys.argv[1:]

for file in files:
    with open(file, 'r') as f:
        jar = f.readlines()

    jar = jar[4:]

    x = []
    a = []; b = []; c = []
    for line in jar:
        tmp = line.split()
        x.append(float(tmp[1]))
        a.append(float(tmp[2]))
        b.append(float(tmp[3]))
        c.append(float(tmp[4]))

    plt.plot(x, a)
    plt.plot(x, b)
    plt.plot(x, c)

plt.show()
