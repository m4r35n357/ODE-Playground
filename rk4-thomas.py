#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, rk4
from math import sin

class Parameters(namedtuple('ParametersType', ['b'])):
    pass

def get_p():
    return Parameters(b=float(argv[8]))

def ode(x, y, z, p):
    return Components(x=sin(y) - p.b * x,
                      y=sin(z) - p.b * y,
                      z=sin(x) - p.b * z)


print(f'RK4: {argv}', file=stderr)
Context.places, skip, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, skip, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
