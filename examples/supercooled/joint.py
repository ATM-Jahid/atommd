#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

xshift = 0
yshift = 1
rdfType = int(sys.argv[1])

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

    x = [xshift + i for i in x]
    ymin = min(y)
    y = [i - ymin + yshift * off for i in y]
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
        y.append(float(line.split()[rdfType]))

    y = [i + yshift * itr for i in y]
    plt.plot(x, y, color=plt.cm.turbo(itr/5), label=f'{temp} (TW)')

files = sys.argv[2:]

for col, off in zip(range(4), [0, 1, 2, 4]):
    csv_draw(files[0], col, off)

for file, itr in zip(files[-1:0:-1], range(5)):
    dot_draw(file, itr)

plt.xlabel('r')
match rdfType:
    case 1:
        plt.ylabel(r'$g_{AA}(r)$')
    case 2:
        plt.ylabel(r'$g_{BB}(r)$')
    case 3:
        plt.ylabel(r'$g_{AB}(r)$')
plt.legend(loc='upper right')
plt.savefig('trial.pdf')
plt.show()
