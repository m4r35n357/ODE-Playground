#!/usr/bin/env python3
#
#  Example: ./tsm-rucklidge.py 15 10 0.01 15000 1 0 0 6.7 2
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, t_sqr, Components


class Parameters(namedtuple('ParametersType', ['alpha', 'kappa'])):
    pass

def get_p(order):
    return Parameters(alpha = float(argv[8]),
                      kappa = float(argv[9]))

def ode(x, y, z, p, i, k):
    return Components(x = p.alpha * y[k] - p.kappa * x[k] - t_prod(y, z, k),
                      y = x[k],
                      z = t_sqr(y, k) - z[k])

tsm(ode, get_p, None)
