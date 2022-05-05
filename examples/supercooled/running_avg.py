#!/usr/bin/env python3

import sys
from statistics import mean, stdev
import matplotlib.pyplot as plt

file = sys.argv[1]

with open(file, 'r') as f:
    jar = f.readlines()

jar.pop()
time = []
tot_e = []; avg_tot_e = []
pot_e = []; avg_pot_e = []
pressure = []; avg_pressure = []

for line in jar:
    tmp = line.split()
    time.append(float(tmp[1]))
    tot_e.append(float(tmp[4]))
    pressure.append(float(tmp[5]))
    pot_e.append(float(tmp[4]) - float(tmp[3]))

for ind in range(len(time)):
    avg_tot_e.append(mean(tot_e[:ind+1]))
    avg_pot_e.append(mean(pot_e[:ind+1]))
    avg_pressure.append(mean(pressure[:ind+1]))

print(f'Potential: {avg_pot_e[-1]}')
print(f'Total: {avg_tot_e[-1]}')
print(f'Pressure: {avg_pressure[-1]}')

plt.plot(time, avg_pot_e, label='Potential Energy')
plt.plot(time, avg_tot_e, label='Total Energy')
plt.plot(time, avg_pressure, label='Pressure')

plt.xlabel('Elapsed time')
plt.ylabel('Dimensionless units')
plt.legend()
plt.show()
