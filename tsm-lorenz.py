#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, Components, Context

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def get_p(order):
    return Parameters(σ=float(argv[8]),
                      ρ=float(argv[9]),
                      β=float(argv[10]) / float(argv[11]))

def ode(x, y, z, p, k):
    return Components(x=p.σ * (y[k] - x[k]),
                      y=p.ρ * x[k] - y[k] - t_prod(x, z, k),
                      z=t_prod(x, y, k) - p.β * z[k])


Context.places, n, δt, n_steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, n, δt, n_steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p(n))
