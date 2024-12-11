#!/usr/bin/env python3
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['a'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]))

def ode(x, y, z, p):
    return Components(x = - p.a * x - 4.0 * (y + z) - y**2,
                      y = - p.a * y - 4.0 * (z + x) - z**2,
                      z = - p.a * z - 4.0 * (x + y) - x**2)

Context.places, skip, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, skip, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
