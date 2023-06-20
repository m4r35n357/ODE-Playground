#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['a'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]))

def ode(x, y, z, p):
    return Components(x=-p.a * x - 4.0 * (y + z) - y * y,
                      y=-p.a * y - 4.0 * (z + x) - z * z,
                      z=-p.a * z - 4.0 * (x + y) - x * x)


Context.places, skip, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, skip, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
