#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, Components

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def ode(x, y, z, p, i, k):
    return Components(x = p.σ * (y[k] - x[k]),
                      y = p.ρ * x[k] - y[k] - t_prod(x, z, k),
                      z = t_prod(x, y, k) - p.β * z[k])

def get_p():
    return Parameters(σ = float(argv[8]),
                      ρ = float(argv[9]),
                      β = float(argv[10]) / float(argv[11]))

tsm(ode, get_p, None)
