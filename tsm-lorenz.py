#!/usr/bin/env python3
#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_mul

def tsm_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def get_p():
    return Parameters(σ=float(argv[8]), ρ=float(argv[9]), β=float(argv[10]) / float(argv[11]))

def ode(x, y, z, p, k):
    return Components(x=p.σ * (y[k] - x[k]),
                      y=p.ρ * x[k] - y[k] - t_mul(x, z, k),
                      z=t_mul(x, y, k) - p.β * z[k])

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
