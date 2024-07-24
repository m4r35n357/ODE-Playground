#!/usr/bin/env python3
# ./e3.py 0.8 0.001 1.5 3 1.701658972 0.760578520 0.430175117
from sys import argv, stderr
from math import log10
import matplotlib.pyplot as plt

def e3(p, w):
    g = 1.0 + 1j * w * p[0]
    g = 1.0 / g + 1.0 / (1.0 / (1j * w * p[1]) + 1j * w * p[2])
    g = 1.0 / g + 1j * w * p[0]
    return 1.0 - abs((g - 1.0) / (g + 1.0))**2

print(f'3rd-order filter frequency response: {argv}', file=stderr)
if len(argv) == 8:
    pb = float(argv[1])
    sb = float(argv[2])
    ksi = float(argv[3])
    dim = int(argv[4])
    values = [0.0] * dim
    for i in range(dim):
        values[i] = float(argv[i + 5])
    L = 1025
    w_axis = [0.0] * L
    response = [0.0] * L
    for i in range(L):
        w_axis[i] = pow(10.0, 2.0 * i / L - 1.0)
        t = e3(values, w_axis[i])
        response[i] = 10.0 * log10(t) if t > 1.0e-18 else -180.0
    plt.figure()
    plt.title(argv, fontsize=10)
    plt.semilogx()
    plt.plot(w_axis, response)
    plt.xlabel('w / w0')
    plt.ylabel('Transmission [dB]')
    y_lim = -80.0
    plt.ylim(y_lim, 0.0)
    plt.hlines(y=10.0 * log10(pb), xmin=1.0e-1, xmax=1.0, colors='r')
    plt.hlines(y=10.0 * log10(sb), xmin=ksi, xmax=1.0e1, colors='r')
    plt.vlines(x=1.0, ymin=y_lim, ymax=0.0, colors='b')
    plt.vlines(x=ksi, ymin=y_lim, ymax=0.0, colors='b')
    plt.grid()
    plt.show()
else:
    raise Exception(f'>>> Wrong number of arguments ({len(argv):8})<<<')

