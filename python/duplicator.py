#!/usr/bin/env python3

import os
import sys

exFile = sys.argv[1]
with open(exFile, 'r') as f:
    jar = f.readlines()

nums = [256, 500, 864, 1372, 2048]
for i in range(2):
    for j in range(2):
        for num in nums:
            jar[2] = str(num) + '\n'
            jar[3] = str(i) + '\n'
            jar[4] = str(j) + '\n'
            writeFile = f'liquid_{num}_{i}_{j}.in'
            with open(writeFile, 'a') as f:
                for line in jar:
                    f.write(line)
            os.system(f'./a.out {writeFile}')
