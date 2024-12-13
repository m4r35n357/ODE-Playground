#!/usr/bin/env python3
#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def get_p():
    return Parameters(σ=float(argv[8]), ρ=float(argv[9]), β=float(argv[10]) / float(argv[11]))

def ode(x, y, z, p):
    return Components(x = p.σ * (y - x),
                      y = (p.ρ - z) * x - y,
                      z = x * y - p.β * z)

Context.places, skip, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, skip, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
