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
        a.append(float(line.split()[1]))
        b.append(float(line.split()[2]))
        c.append(float(line.split()[3]))

    d = [0.2*i + 0.8*j for i, j in zip(a, b)]
    plt.plot(x, c, label=f'{temp} ab')
    plt.plot(x, d, label=f'{temp} fo')

files = sys.argv[1:]

for file, itr in zip(files, range(5)):
    print(file, itr)
    dot_draw(file, itr)

plt.legend()
plt.show()
