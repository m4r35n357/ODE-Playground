#!/usr/bin/env python3

from sys import argv
from taylor import jet_0, t_prod, t_horner, jet_c, t_sqr

n, h, steps = int(argv[1]), float(argv[2]), int(argv[3])
x, y, z = float(argv[4]), float(argv[5]), float(argv[6])

cx, cy, cz = jet_0(n + 1), jet_0(n + 1), jet_0(n + 1)
w1 = jet_c(1.0, n)
for step in range(1, steps):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
    cx[0], cy[0], cz[0] = x, y, z
    for k in range(n):
        cx[k + 1] = (cy[k] + 2.0 * t_prod(cx, cy, k) + t_prod(cx, cz, k)) / (k + 1)
        cy[k + 1] = (w1[k] - 2.0 * t_sqr(cx, k) + t_prod(cy, cz, k)) / (k + 1)
        cz[k + 1] = (cx[k] - t_sqr(cx, k) - t_sqr(cy, k)) / (k + 1)
    x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
