#!/usr/bin/env python3
#
#  Example: ./rk4-lorenz.py  6 1  .01 10000  -15.8 -17.48 35.64  10 28 8 3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def get_p():
    return Parameters(σ=float(argv[8]), ρ=float(argv[9]), β=float(argv[10]) / float(argv[11]))

def ode(x, y, z, p):
    return Components(x=p.σ * (y - x),
                      y=p.ρ * x - y - x * z,
                      z=x * y - p.β * z)


Context.places, n, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, n, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
