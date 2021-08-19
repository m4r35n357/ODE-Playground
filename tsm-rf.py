#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_prod, Components, t_jet, t_sqr

class Parameters(namedtuple('ParametersType', ['α', 'γ', 'jet1', 'a', 'b', 'c'])):
    pass

class Intermediates(namedtuple('IntermediatesType', [])):
    pass

def get_p(order):
    return Parameters(α=t_jet(order, float(argv[9])), γ=float(argv[10]),
                      jet1=t_jet(order, 1.0), a=t_jet(order), b=t_jet(order), c=t_jet(order))

def ode(x, y, z, p, k):
    p.a[k] = z[k] + t_sqr(x, k) - p.jet1[k]
    p.b[k] = 4.0 * z[k] - p.a[k]
    p.c[k] = p.α[k] + t_prod(x, y, k)
    return Components(x=t_prod(y, p.a, k) + p.γ * x[k],
                      y=t_prod(x, p.b, k) + p.γ * y[k],
                      z=- 2.0 * t_prod(z, p.c, k))


tsm(ode, get_p)
