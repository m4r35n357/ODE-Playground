from sys import argv
from math import log10
import numpy as np
import matplotlib.pyplot as plt

def tx(n, p, w):
    g = complex(1.0, w * p[0])
    for c in range(1, n):
        g = 1.0 / g + complex(0.0, w * p[c])
    for c in range(n - 2, -1, -1):
        g = 1.0 / g + complex(0.0, w * p[c])
    return 1.0 - abs((g - 1.0) / (g + 1.0))**2

dim = int(argv[1])
values = np.zeros(dim)
for i in range(dim):
    values[i] = float(argv[i + 2])
L = 1024
w_axis = np.logspace(-3, 3, num=L)
response = np.zeros(L)
for i in range(L):
    t = tx(dim, values, w_axis[i])
    response[i] = 10.0 * log10(t) if t > 1.0e-18 else -180.0

plt.figure()
plt.semilogx()
#plt.plot(np.linspace(0, 10, len(response)), response)
plt.plot(w_axis, response)
plt.xlabel('Normalized frequency')
plt.ylabel('Gain [dB]')
plt.grid()
plt.show()