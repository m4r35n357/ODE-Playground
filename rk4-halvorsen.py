#!/usr/bin/env python3
#
#  Example: ./rk4-halvorsen.py  6 1  .01 10000  1 0 0  1.4 4
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['a', 'b'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]), b=float(argv[9]))

def ode(x, y, z, p):
    return Components(x=-p.a * x - p.b * (y + z) - y * y,
                      y=-p.a * y - p.b * (z + x) - z * z,
                      z=-p.a * z - p.b * (x + y) - x * x)


Context.places, n, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, n, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
