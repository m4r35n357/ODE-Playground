#!/usr/bin/env python3
# ./bw.py 0.8 0.0001 1.2 6 1.569771623 0.931740540 0.344074589 2.094462346 0.532291777 1.233067458
from sys import argv, stderr
from math import log10
import matplotlib.pyplot as plt

def bw(n, p, w):
    g = complex(1.0, w * p[0])
    for r in range(n):
        g = 1.0 / g + complex(0.0, w * p[r])
    for r in range(n - 2, -1, -1):
        g = 1.0 / g + complex(0.0, w * p[r])
    return 1.0 - abs((g - 1.0) / (g + 1.0))**2

print(f'Butterworth filter frequency response: {argv}', file=stderr)
if len(argv) == 8:
    dim = int(argv[1])
    values = [0.0] * dim
    for i in range(dim):
        values[i] = float(argv[i + 2])
    L = 1025
    w_axis = [0.0] * L
    response = [0.0] * L
    for i in range(L):
        w_axis[i] = pow(10.0, 2.0 * i / L - 1.0)
        t = bw(dim, values, w_axis[i])
        response[i] = 10.0 * log10(t) if t > 1.0e-18 else -180.0
    plt.figure()
    plt.title(argv, fontsize=10)
    plt.semilogx()
    plt.plot(w_axis, response)
    plt.xlabel('w / w0')
    plt.ylabel('Transmission [dB]')
    y_lim = -80.0
    plt.ylim(y_lim, 0.0)
    plt.grid()
    plt.show()
else:
    raise Exception(f'>>> Wrong number of arguments ({len(argv):11})<<<')
