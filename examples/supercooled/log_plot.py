#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

def dot_draw(fileName, itr):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    jar = jar[2:]

    x = []; y = []
    for line in jar:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[8]))

    x = [(i-50000)*0.01 for i in x]
    plt.plot(x, y)

files = sys.argv[1:]

for file, itr in zip(files, range(5)):
    print(file, itr)
    dot_draw(file, itr)

plt.xscale('log')
plt.yscale('log')
plt.show()
