#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from collections import namedtuple
from ad import tsm, t_sin_cos, Components, t_jet

class Parameters(namedtuple('ParametersType', ['b', 'sx', 'sy', 'sz', 'cx', 'cy', 'cz'])):
    pass

def get_p(order):
    return Parameters(b=float(argv[8]),
                      sx=t_jet(order), cx=t_jet(order),
                      sy=t_jet(order), cy=t_jet(order),
                      sz=t_jet(order), cz=t_jet(order))

def ode(x, y, z, p, k):
    return Components(x=t_sin_cos(p.sy, p.cy, y, k)[0] - p.b * x[k],
                      y=t_sin_cos(p.sz, p.cz, z, k)[0] - p.b * y[k],
                      z=t_sin_cos(p.sx, p.cx, x, k)[0] - p.b * z[k])


tsm(ode, get_p)
