#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, Components, t_jet, t_sqr

class Parameters(namedtuple('ParametersType', ['α', 'γ'])):
    pass

class Intermediates(namedtuple('IntermediatesType', ['jet1', 'a', 'b', 'c'])):
    pass

def ode(x, y, z, p, i, k):
    i.a[k] = z[k] + t_sqr(x, k) - i.jet1[k]
    i.b[k] = 4.0 * z[k] - i.a[k]
    i.c[k] = p.α[k] + t_prod(x, y, k)
    return Components(x = t_prod(y, i.a, k) + p.γ * x[k],
                      y = t_prod(x, i.b, k) + p.γ * y[k],
                      z = - 2.0 * t_prod(z, i.c, k))

def get_p(order):
    return Parameters(α = t_jet(order, float(argv[9])),
                      γ = float(argv[10]))

def get_i(order):
    return Intermediates(jet1 = t_jet(order, 1.0),
                         a = t_jet(order),
                         b = t_jet(order),
                         c = t_jet(order))

tsm(ode, get_p, get_i)
