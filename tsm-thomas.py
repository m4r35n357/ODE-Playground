#!/usr/bin/env python3
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, tsm_jet, t_sin_cos

class Parameters(namedtuple('ParametersType', ['b', 'sx', 'sy', 'sz', 'cx', 'cy', 'cz'])):
    pass

def get_p(n):
    return Parameters(b=float(argv[8]), sx=tsm_jet(n), cx=tsm_jet(n), sy=tsm_jet(n), cy=tsm_jet(n), sz=tsm_jet(n), cz=tsm_jet(n))

def ode(x, y, z, p, k):
    p.sx[k], p.cx[k] = t_sin_cos(p.sx, p.cx, x, k)
    p.sy[k], p.cy[k] = t_sin_cos(p.sy, p.cy, y, k)
    p.sz[k], p.cz[k] = t_sin_cos(p.sz, p.cz, z, k)
    return Components(x=p.sy[k] - p.b * x[k],
                      y=p.sz[k] - p.b * y[k],
                      z=p.sx[k] - p.b * z[k])

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p(order))
