#!/usr/bin/env python3

from sys import argv
from taylor import jet_0, t_horner, t_sqr

n, h, steps = int(argv[1]), float(argv[2]), int(argv[3])
x, y, z = float(argv[4]), float(argv[5]), float(argv[6])
a = float(argv[7])

cx, cy, cz = jet_0(n + 1), jet_0(n + 1), jet_0(n + 1)
for step in range(1, steps):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
    cx[0], cy[0], cz[0] = x, y, z
    for k in range(n):
        cx[k + 1] = - (a * cx[k] + 4.0 * cy[k] + 4.0 * cz[k] + t_sqr(cy, k)) / (k + 1)
        cy[k + 1] = - (a * cy[k] + 4.0 * cz[k] + 4.0 * cx[k] + t_sqr(cz, k)) / (k + 1)
        cz[k + 1] = - (a * cz[k] + 4.0 * cx[k] + 4.0 * cy[k] + t_sqr(cx, k)) / (k + 1)
    x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
