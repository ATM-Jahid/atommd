#!/usr/bin/env python3

import sys
from statistics import mean
import matplotlib.pyplot as plt

column = int(sys.argv[1])

def dot_draw(fileName):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    jar = jar[-100:]

    x = []; y = []
    for line in jar:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[column]))

    m = mean(y)
    print(m)
    plt.plot(x, y)
    plt.plot(x, [m for i in x])

files = sys.argv[2:]

for file in files:
    print(file)
    dot_draw(file)

plt.show()
