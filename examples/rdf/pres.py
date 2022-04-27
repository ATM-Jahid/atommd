#!/usr/bin/env python3

import sys
from statistics import mean

file = sys.argv[1]

with open(file, 'r') as f:
    jar = f.readlines()
jar.pop()

pres = []
for line in jar:
    pres.append(float(line.split()[5]))

print(mean(pres))
print(mean(pres[200:]))
