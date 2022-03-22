#!/usr/bin/env python3

import matplotlib.pyplot as plt

with open('walltimes.csv') as f:
    jar = f.readlines()

jar.pop(0)
num = []
allp = []
cell = []
neigh = []
nebr_cell= []

for line in jar:
    tmp = line.split(',')
    num.append(int(tmp[0]))
    allp.append(int(tmp[1]))
    cell.append(int(tmp[2]))
    neigh.append(int(tmp[3]))
    nebr_cell.append(int(tmp[4]))

plt.plot(num, allp, label='all_pair')
plt.plot(num, cell, label='cell_sub')
plt.plot(num, neigh, label='neigh')
plt.plot(num, nebr_cell, label='neigh_cell')

plt.xlabel('Number of atoms')
plt.ylabel('Seconds')

plt.legend()
plt.savefig('wall.pdf')
