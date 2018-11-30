#!/usr/bin/env python3

from sys import argv
from taylor import jet_0, t_horner, t_sin_cos

n, h, steps = int(argv[1]), float(argv[2]), int(argv[3])
x, y, z = float(argv[4]), float(argv[5]), float(argv[6])
b = float(argv[7])

cx, cy, cz = jet_0(n + 1), jet_0(n + 1), jet_0(n + 1)
wsx, wcx = jet_0(n), jet_0(n)
wsy, wcy = jet_0(n), jet_0(n)
wsz, wcz = jet_0(n), jet_0(n)
for step in range(1, steps):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
    cx[0], cy[0], cz[0] = x, y, z
    for k in range(n):
        wsx[k], wcx[k] = t_sin_cos(wsx, wcx, cx, k)
        wsy[k], wcy[k] = t_sin_cos(wsy, wcy, cy, k)
        wsz[k], wcz[k] = t_sin_cos(wsz, wcz, cz, k)
        cx[k + 1] = (wsy[k] - b * cx[k]) / (k + 1)
        cy[k + 1] = (wsz[k] - b * cy[k]) / (k + 1)
        cz[k + 1] = (wsx[k] - b * cz[k]) / (k + 1)
    x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
