#!/usr/bin/env python3
#
#  Example: ./rk4-lorenz.py 12 NA 1 .01 100000 -15.8 -17.48 35.64 10 28 8 3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import rk4, t_prod, Components

class Parameters(namedtuple('ParametersType', ['σ', 'ρ', 'β'])):
    pass

def get_p():
    return Parameters(σ = float(argv[9]),
                      ρ = float(argv[10]),
                      β = float(argv[11]) / float(argv[12]))

def ode(x, y, z, p):
    return Components(x = p.σ * (y - x),
                      y = p.ρ * x - y - x * z,
                      z = x * y - p.β * z)

rk4(ode, get_p)
