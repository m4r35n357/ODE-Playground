#!/usr/bin/env python3
#
#  Example: ./tsm-halvorsen.py  6 8  .01 10000  1 0 0  1.4 4
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, tsm, t_sqr

class Parameters(namedtuple('ParametersType', ['a', 'b'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]), b=float(argv[9]))

def ode(x, y, z, p, k):
    return Components(x=-p.a * x[k] - p.b * (y[k] + z[k]) - t_sqr(y, k),
                      y=-p.a * y[k] - p.b * (z[k] + x[k]) - t_sqr(z, k),
                      z=-p.a * z[k] - p.b * (x[k] + y[k]) - t_sqr(x, k))


Context.places, n, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, n, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())