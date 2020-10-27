#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Example: ./tsm-thomas.py 15 10 0.1 30000 1 0 0 .19
#

from sys import argv
from collections import namedtuple
from ad import tsm, Components, t_jet, t_sin_cos

class Parameters(namedtuple('ParametersType', ['b'])):
    pass

class Intermediates(namedtuple('IntermediatesType', ['sx', 'cx', 'sy', 'cy', 'sz', 'cz'])):
    pass

def get_p(order):
    return Parameters(b = float(argv[8]))

def get_i(order):
    return Intermediates(sx = t_jet(order), cx = t_jet(order),
                         sy = t_jet(order), cy = t_jet(order),
                         sz = t_jet(order), cz = t_jet(order))

def ode(x, y, z, p, i, k):
    i.sx[k], i.cx[k] = t_sin_cos(i.sx, i.cx, x, k)
    i.sy[k], i.cy[k] = t_sin_cos(i.sy, i.cy, y, k)
    i.sz[k], i.cz[k] = t_sin_cos(i.sz, i.cz, z, k)
    return Components(x = i.sy[k] - p.b * x[k],
                      y = i.sz[k] - p.b * y[k],
                      z = i.sx[k] - p.b * z[k])

tsm(ode, get_p, get_i)
