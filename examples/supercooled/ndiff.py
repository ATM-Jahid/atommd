#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

def dot_draw(fileName, itr):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    jar = jar[2:]

    for jj in [8,9,10]:
        x = []; y = []
        for line in jar:
            x.append(float(line.split()[0]))
            y.append(float(line.split()[jj]))

        x = [(i-50000)*0.001 for i in x]
        m, b = np.polyfit(x, y, 1)
        coords = [m*i+b for i in x]
        print(m/6)
        plt.plot(x, y)
        plt.plot(x, coords)

files = sys.argv[1:]

for file, itr in zip(files, range(5)):
    print(file, itr)
    dot_draw(file, itr)

plt.show()
