#!/usr/bin/env python3

import sys
import itertools
import matplotlib.pyplot as plt

def draw(fileName, mark, style, hello):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    jar.pop(0)
    temp = []
    pot = []
    tot = []
    pres = []
    for line in jar:
        tmp = line.split(',')
        temp.append(float(tmp[0]))
        pot.append(float(tmp[1]))
        tot.append(float(tmp[2]))
        pres.append(float(tmp[3]))

    overT = [1/x for x in temp]
    pres = [x/10 for x in pres]

    plt.plot(overT, pot, color=plt.cm.turbo(0.9),
            marker=next(mark), ls=style, label=f'Potential E. ({hello})')
    plt.plot(overT, tot, color=plt.cm.turbo(0.6),
            marker=next(mark), ls=style, label=f'Total E. ({hello})')
    plt.plot(overT, pres, color=plt.cm.turbo(0.3),
            marker=next(mark), ls=style, label=f'Pressure ({hello})')

files = sys.argv[1:]
markers1 = itertools.cycle(('o', 's', 'p'))
markers2 = itertools.cycle(('v', '^', '>'))
draw(files[0], markers1, '-', 'this work')
draw(files[1], markers2, ':', 'Kob-Andersen')

plt.xlabel(r'1/T')
plt.ylabel(r'$E_{pot}$, $E_{tot}$, $p$')
plt.legend()
plt.savefig('enerKob.pdf')
plt.show()
