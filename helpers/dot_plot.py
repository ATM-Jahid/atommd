#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

files = sys.argv[1:]

for file in files:
    with open(file, 'r') as f:
        jar = f.readlines()

    jar = jar[1:]

    x = []; y = []
    for line in jar:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

    plt.plot(x, y)
    plt.show()
