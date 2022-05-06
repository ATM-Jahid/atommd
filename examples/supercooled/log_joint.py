#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

msdType = int(sys.argv[1])

def csv_draw(fileName, column, off):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    temp = float(jar[0].split(',')[column+1])
    jar = jar[1:]

    x = []; y = []
    for line in jar:
        tmp = line.split(',')
        x.append(float(tmp[0]))
        y.append(float(tmp[column+1]))

    x = [i / 48**0.5 for i in x]
    plt.plot(x, y, color=plt.cm.turbo(off/5), ls=':', label=f'{temp} (KA)')

def dot_draw(fileName, itr):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    ind = fileName.find('_') + 1
    temp = float(fileName[ind:ind+3])
    jar = jar[2:]

    x = []; y = []
    for line in jar:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[7+msdType]))

    x = [(i-50000)*0.001 for i in x]
    plt.plot(x, y, color=plt.cm.turbo(itr/5), label=f'{temp} (TW)')

files = sys.argv[2:]

#for col, off in zip(range(4), [0, 1, 2, 4]):
#    csv_draw(files[0], col, off)

for file, itr in zip(files[-1:0:-1], range(5)):
    dot_draw(file, itr)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('time')
plt.ylabel(r'$\langle r^2 (t) \rangle_{BB}$')
plt.legend()
plt.savefig('trial.pdf')
plt.show()
