#!/usr/bin/env python3
#
#  Example: ./tsm-halvorsen.py 15 _ 10 .01 100001 1 0 0 1.4
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, Components, t_sqr

class Parameters(namedtuple('ParametersType', ['α'])):
    pass

def get_p(order):
    α = float(argv[9])
    return Parameters(α = float(argv[9]))

def ode(x, y, z, p, i, k):
    return Components(x = - p.α * x[k] - 4.0 * (y[k] + z[k]) - t_sqr(y, k),
                      y = - p.α * y[k] - 4.0 * (z[k] + x[k]) - t_sqr(z, k),
                      z = - p.α * z[k] - 4.0 * (x[k] + y[k]) - t_sqr(x, k))

tsm(ode, get_p, None)
