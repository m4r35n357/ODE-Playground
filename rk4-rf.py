#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['α', 'γ'])):
    pass

def get_p(order):
    return Parameters(α=float(argv[8]), γ=float(argv[9]))

def ode(x, y, z, p):
    a = z + x * x - 1.0
    b = 4.0 * z - a
    c = p.α + x * y
    return Components(x=y * a + p.γ * x,
                      y=x * b + p.γ * y,
                      z=- 2.0 * z * c)


Context.places, n, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, n, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p(n))
