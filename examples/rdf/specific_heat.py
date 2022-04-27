#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

m = 2.8258994058129354
z = -8.235216022681811
a = 1.8705138521250742
b = -7.331111553575411
c = 6.6610516401005775

lin_x = np.linspace(0.2, 2.4, 100)
lin_y = [m for x in lin_x]
plt.plot(lin_x, lin_y,
        label='Before melting', color=plt.cm.YlOrRd(0.4))

quad_x = np.linspace(2.6, 3.0, 20)
quad_y = [2*a*x+b for x in quad_x]
plt.plot(quad_x, quad_y,
        label='After melting', color=plt.cm.YlOrRd(0.8))

plt.xlabel('Temperature')
plt.ylabel('Specific heat')
plt.legend()
plt.savefig('specHeat.pdf')
plt.show()
