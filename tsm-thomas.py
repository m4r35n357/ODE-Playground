#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_jet, t_sin_cos

class Parameters(namedtuple('ParametersType', ['b', 'sx', 'sy', 'sz', 'cx', 'cy', 'cz'])):
    pass

def get_p(n):
    return Parameters(b=float(argv[8]), sx=t_jet(n), cx=t_jet(n), sy=t_jet(n), cy=t_jet(n), sz=t_jet(n), cz=t_jet(n))

def ode(x, y, z, p, k):
    return Components(x=t_sin_cos(p.sy, p.cy, y, k)[0] - p.b * x[k],
                      y=t_sin_cos(p.sz, p.cz, z, k)[0] - p.b * y[k],
                      z=t_sin_cos(p.sx, p.cx, x, k)[0] - p.b * z[k])


print(f'TSM: {argv}', file=stderr)
Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p(order))
