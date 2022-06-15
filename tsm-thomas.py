#!/usr/bin/env python3
#
#  Example: ./tsm-thomas.py 6 8 0.1 30000 1 0 0 .19
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, tsm, t_jet, t_sin_cos

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


Context.places, n, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, n, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p(n))
