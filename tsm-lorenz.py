#!/usr/bin/env python3
#
#  Example: ./tsm-lorenz.py 6 8 .01 10000 -15.8 -17.48 35.64 10 28 8 3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, tsm, t_prod

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def get_p():
    return Parameters(σ=float(argv[8]), ρ=float(argv[9]), β=float(argv[10]) / float(argv[11]))

def ode(x, y, z, p, k):
    return Components(x=p.σ * (y[k] - x[k]),
                      y=p.ρ * x[k] - y[k] - t_prod(x, z, k),
                      z=t_prod(x, y, k) - p.β * z[k])


Context.places, n, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, n, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
