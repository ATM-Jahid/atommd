#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

folder = sys.argv[1:]

def draw(files, wn, sv, pv, tag, mark):
    temps = []
    diffs = []
    for file in files:
        temp = float(file[3:6])
        temps.append(temp)
        with open(file, 'r') as f:
            jar = f.readlines()
        diffs.append(float(jar[-1].split()[wn]))

    overT = [1/T for T in temps]
    k = 5.732e-8
    diffuse = [k*x for x in diffs]
    plt.scatter(overT, diffuse,
            color=plt.cm.turbo(sv), marker=mark, label=tag)

    m, b = np.polyfit(overT, np.log(diffuse), 1)
    ordinates = [m*x+b for x in overT]
    plt.plot(overT, np.exp(ordinates), color=plt.cm.turbo(pv))

msdFiles = [x for x in folder if 'dfs' in x]
vacFiles = [x for x in folder if 'acf' in x]

draw(msdFiles, 1, 0.9, 0.9, 'From MSD', 'o')
draw(vacFiles, 2, 0.2, 0.2, 'From VACF', 'v')

plt.yscale('log')
plt.xlabel('1/T')
plt.ylabel('Diffusion coefficient')
plt.tight_layout()
plt.legend()
plt.savefig('diffComp.pdf')
plt.show()
