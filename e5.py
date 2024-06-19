from sys import argv
from math import log10
import matplotlib.pyplot as plt

def e5(n, p, w):
    g = complex(1.0, w * p[0])
    g = 1.0 / g + 1.0 / (complex(0.0, w * p[2]) + 1.0 / complex(0.0, w * p[1]))
    g = 1.0 / g + complex(0.0, w * p[3])
    g = 1.0 / g + 1.0 / (complex(0.0, w * p[5]) + 1.0 / complex(0.0, w * p[4]))
    g = 1.0 / g + complex(0.0, w * p[0])
    return 1.0 - abs((g - 1.0) / (g + 1.0))**2

dim = int(argv[1])
values = [0.0] * dim
for i in range(dim):
    values[i] = float(argv[i + 2])
L = 1025
w_axis = [0.0] * L
response = [0.0] * L
for i in range(L):
    w_axis[i] = pow(10.0, 5.0 * (2.0 * i / L - 1.0))
    t = e5(dim, values, w_axis[i])
    response[i] = 10.0 * log10(t) if t > 1.0e-18 else -180.0

if len(argv) == 8:
    plt.figure()
    plt.semilogx()
    plt.plot(w_axis, response)
    plt.xlabel('w / w0')
    plt.ylabel('Transmission [dB]')
    plt.grid()
    plt.show()
else:
    raise Exception(f'>>> Wrong number of arguments ({len(argv):8})<<<')
