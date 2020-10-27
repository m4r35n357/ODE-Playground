#!/usr/bin/env python3
#
#  Example: ./tsm-rossler.py 15 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, Components, t_jet

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c'])):
    pass

def get_p(order):
    return Parameters(a = float(argv[8]),
                      b = t_jet(order, float(argv[9])),
                      c = float(argv[10]))

def ode(x, y, z, p, i, k):
    return Components(x = - y[k] - z[k],
                      y = x[k] + p.a * y[k],
                      z = p.b[k] + t_prod(x, z, k) - p.c * z[k])

tsm(ode, get_p, None)
