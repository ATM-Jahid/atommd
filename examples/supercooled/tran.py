#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

files = sys.argv[1:]

for file in files:
    with open(file, 'r') as f:
        jar = f.readlines()

    jar.pop(0)
    overT = []
    pot = []
    tot = []
    pres = []
    for line in jar:
        tmp = line.split(',')
        overT.append(float(tmp[0]))
        pot.append(float(tmp[1]))
        tot.append(float(tmp[2]))
        pres.append(float(tmp[3]))

    temp = [1/x for x in overT]
    pres = [x*10 for x in pres]

    with open('fig1tran.csv', 'w') as f:
        f.write('temp,potential,total,pressure\n')
        for x, p, t, pr in zip(temp, pot, tot, pres):
            f.write(f'{x},{p},{t},{pr}\n')
