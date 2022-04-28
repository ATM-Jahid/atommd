#!/usr/bin/env python3

import os
import sys

exFile = sys.argv[1]
with open(exFile, 'r') as f:
    jar = f.readlines()

temps = [x/10 for x in range(7, 15)]
for temp in temps:
    jar[0] = str(temp) + '\n'
    writeFile = f'Ar_{temp}.in'
    with open(writeFile, 'w') as f:
        for line in jar:
            f.write(line)
    os.system(f'./a.out {writeFile}')
