#!/usr/bin/env python3

import sys
import numpy as np
from statistics import mean, stdev
import matplotlib.pyplot as plt

files = sys.argv[1:]

temps = []
avg_tot = []
avg_kin = []

for file in files:
    temp = float(file[3:6])
    temps.append(temp)

    with open(file, 'r') as f:
        jar = f.readlines()

    jar.pop()
    tot_e = []
    kin_e = []

    for line in jar:
        tmp = line.split()
        kin_e.append(float(tmp[3]))
        tot_e.append(float(tmp[4]))

    avg_tot.append(mean(tot_e))
    avg_kin.append(mean(kin_e))

#plt.scatter(temps, avg_kin, label='Kinetic Energy')
#plt.scatter(temps, avg_tot, label='Total Energy')

lin_energy = avg_tot[:-3]
plt.scatter(temps[:-3], lin_energy, s=15,
        label='Before melting', color=plt.cm.YlOrRd(0.4))
ml, bl = np.polyfit(temps[:-3], lin_energy, 1)
print(ml, bl)
ordinates_lin = [ml*x+bl for x in temps[:-3]]
plt.plot(temps[:-3], ordinates_lin, color=plt.cm.YlOrRd(0.4))

quad_energy = avg_tot[-3:]
plt.scatter(temps[-3:], quad_energy, marker='v', s=15,
        label='After melting', color=plt.cm.YlOrRd(0.8))
aq, bq, cq = np.polyfit(temps[-3:], quad_energy, 2)
print(aq, bq, cq)
fine_temp = np.linspace(temps[-3], temps[-1], 20)
ordinates_quad = [aq*x**2+bq*x+cq for x in fine_temp]
plt.plot(fine_temp, ordinates_quad, color=plt.cm.YlOrRd(0.8))

plt.xlabel('Temperature')
plt.ylabel('Energy')
plt.legend()
plt.savefig('EvsT.pdf')
plt.show()
