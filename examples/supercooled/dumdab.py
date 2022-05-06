#!/usr/bin/env python3

import sys

files = sys.argv[1:]

def extract(fileName):
    with open(fileName, 'r') as f:
        jar = f.readlines()

    ini = jar[50*873+9:51*873]
    end = jar[149*873+9:150*873]

    xi = 0; xf = 0
    yi = 0; yf = 0
    zi = 0; zf = 0

    for line in ini:
        tmp = line.split()
        if tmp[1] == '1':
            xi += float(tmp[2])
            yi += float(tmp[3])
            zi += float(tmp[4])

    for line in end:
        tmp = line.split()
        if tmp[1] == '1':
            xf += float(tmp[2])
            yf += float(tmp[3])
            zf += float(tmp[4])

    mR = 1
    Q = 1/0.8/0.2 * (mR*0.8+0.2)**2
    D = Q * ((xf-xi)**2 + (yf-yi)**2 + (zf-zi)**2) / 6 / 864 / 100

    print(fileName)
    print(D)

for file in files:
    extract(file)
