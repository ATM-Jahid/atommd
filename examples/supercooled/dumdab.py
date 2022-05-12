#!/usr/bin/env python3

import sys
from statistics import mean
import matplotlib.pyplot as plt

files = sys.argv[1:]

def extract(fileName, tint):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    ind = fileName.find('_')
    temp = fileName[ind+1:ind+4]

    buffDiff = 5
    buffLen = 50
    buff = [0] * buffLen

    for b in range(0, 101-buffLen, buffDiff):
        dabs = []
        ini = jar[(50+b)*873+9:(51+b)*873]

        for jj in range(1, 51):
            end = jar[(50+b+jj)*873+9:(50+b+jj+1)*873]

            xi = 0; xf = 0
            yi = 0; yf = 0
            zi = 0; zf = 0

            for line in ini:
                tmp = line.split()
                if tmp[1] == '2':
                    xi += float(tmp[2])
                    yi += float(tmp[3])
                    zi += float(tmp[4])

            for line in end:
                tmp = line.split()
                if tmp[1] == '2':
                    xf += float(tmp[2])
                    yf += float(tmp[3])
                    zf += float(tmp[4])

            mR = 1
            num_atoms = 864
            Q = 1/0.8/0.2 * (mR*0.8+0.2)**2
            f = Q / 6 / num_atoms / jj
            D = f * ((xf-xi)**2 + (yf-yi)**2 + (zf-zi)**2)
            dabs.append(D)

        for ii in range(len(dabs)):
            buff[ii] += dabs[ii]

    numBuff = (100-buffLen) // buffDiff + 1
    buff = [x/numBuff for x in buff]
    runner = [mean(buff[:i+1]) for i in range(50)]
    plt.plot(runner, color=plt.cm.turbo(tint), label=temp)

tints = [i/10 for i in range(1, 10, 2)]
for file, tint in zip(files, tints):
    extract(file, tint)

plt.legend()
plt.show()
