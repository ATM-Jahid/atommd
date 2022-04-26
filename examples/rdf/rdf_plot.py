#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

files = sys.argv[1:]

for file in files:
    temp = float(file[3:6])

    with open(file, 'r') as f:
        jar = f.readlines()

    jar = jar[1:]

    x = []; y = []
    for line in jar:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

    tint = plt.cm.YlOrRd((temp-0.2)/2.8*0.8+0.1)
    #tint = plt.cm.YlOrRd((temp-2.2)*2)
    plt.plot(x, y, label=f'{temp}', color=tint)

plt.xlabel('Radial distance')
plt.ylabel('RDF')
plt.legend()
plt.savefig('ArRdf.pdf')
#plt.savefig('MeltRdf.pdf')
plt.show()
