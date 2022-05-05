#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

def csv_draw(fileName, column):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    temp = float(jar[0].split(',')[column+1])
    jar = jar[1:]

    x = []; y = []
    for line in jar:
        tmp = line.split(',')
        x.append(float(tmp[0]))
        y.append(float(tmp[column+1]))

    plt.plot(x, y, label=f'{temp}')

file = sys.argv[1]

for col in range(4):
    csv_draw(file, col)

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
