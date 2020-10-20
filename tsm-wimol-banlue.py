#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Example: ./tsm-wimol-banlue.py 15 _ 10 0.1 10001 1 0 0 2.0
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, Components, t_jet, t_tan_sec2, t_abs


class Parameters(namedtuple('ParametersType', ['α'])):
    pass

class Intermediates(namedtuple('IntermediatesType', ['tx', 'sx'])):
    pass

def get_p(order):
    return Parameters(α = t_jet(order, float(argv[9])))

def get_i(order):
    return Intermediates(tx = t_jet(order),
                         sx = t_jet(order))

def ode(x, y, z, p, i, k):
    i.tx[k], i.sx[k] = t_tan_sec2(i.tx, i.sx, x, k, hyp=True)
    return Components(x = y[k] - x[k],
                      y = - t_prod(z, i.tx, k),
                      z = - p.α[k] + t_prod(x, y, k) + t_abs(y, k))

tsm(ode, get_p, get_i)
