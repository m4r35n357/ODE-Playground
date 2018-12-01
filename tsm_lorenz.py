#!/usr/bin/env python3

#  Example: ./tsm_lorenz.py 16 10 .01 100001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py

from sys import argv
from taylor import jet_0, t_prod, t_horner

n, h, steps = int(argv[2]), float(argv[3]), int(argv[4])
x, y, z = float(argv[5]), float(argv[6]), float(argv[7])
s, r, b = float(argv[8]), float(argv[9]), float(argv[10]) / float(argv[11])

cx, cy, cz = jet_0(n + 1), jet_0(n + 1), jet_0(n + 1)
for step in range(1, steps):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
    cx[0], cy[0], cz[0] = x, y, z
    for k in range(n):
        cx[k + 1] = s * (cy[k] - cx[k]) / (k + 1)
        cy[k + 1] = (r * cx[k] - t_prod(cx, cz, k) - cy[k]) / (k + 1)
        cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
    x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
