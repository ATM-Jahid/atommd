#!/usr/bin/env python3

import sys
from statistics import mean, stdev
import matplotlib.pyplot as plt
import itertools

files = sys.argv[1:]
print(files)

for file in files:
    with open(file, 'r') as f:
        jar = f.readlines()

    jar.pop()
    time = []
    tot_e = []; avg_tot_e = []; std_tot_e = []
    temperature = []; avg_temperature = []; std_temperature = []
    pressure = []; avg_pressure = []; std_pressure = []

    for line in jar:
        tmp = line.split()
        time.append(float(tmp[1]))
        tot_e.append(float(tmp[4]))
        temperature.append(2*float(tmp[3])/3)
        pressure.append(float(tmp[5]))

    avg_tot_e.append(tot_e[0])
    std_tot_e.append(0)
    avg_temperature.append(temperature[0])
    std_temperature.append(0)
    avg_pressure.append(pressure[0])
    std_pressure.append(0)
    for ind in range(1, len(time)):
        avg_tot_e.append(mean(tot_e[:ind+1]))
        std_tot_e.append(stdev(tot_e[:ind+1]))
        avg_temperature.append(mean(temperature[:ind+1]))
        std_temperature.append(stdev(temperature[:ind+1]))
        avg_pressure.append(mean(pressure[:ind+1]))
        std_pressure.append(stdev(pressure[:ind+1]))

    tag = file[file.find('_')+1:file.find('.')]
    #plt.plot(time, avg_tot_e, label=tag)
    #plt.plot(time[1:], std_tot_e[1:], label=tag)
    #plt.plot(time, avg_temperature, label=tag)
    #plt.plot(time[1:], std_temperature[1:], label=tag)
    #plt.plot(time, avg_pressure, label=tag)
    plt.plot(time[1:], std_pressure[1:], label=tag)

plt.xlabel('Elapsed time')
#plt.ylabel('Dimensionless units')
plt.ylabel('Standard deviation')
#plt.legend(loc=(0.6,0.33))
plt.legend()
plt.savefig('pressure_std.pdf')
#plt.show()
