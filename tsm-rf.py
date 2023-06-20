#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from sys import argv
from collections import namedtuple
from ad import Components, Context, tsm, t_jet, t_const, t_prod, t_sqr

class Parameters(namedtuple('ParametersType', ['α', 'γ', 'a', 'b', 'c'])):
    pass

def get_p(n):
    return Parameters(α=float(argv[8]), γ=float(argv[9]), a=t_jet(n), b=t_jet(n), c=t_jet(n))

def ode(x, y, z, p, k):
    p.a[k] = z[k] + t_sqr(x, k) - t_const(1.0, k)
    p.b[k] = 4.0 * z[k] - p.a[k]
    p.c[k] = t_const(p.α, k) + t_prod(x, y, k)
    return Components(x=t_prod(y, p.a, k) + p.γ * x[k],
                      y=t_prod(x, p.b, k) + p.γ * y[k],
                      z=- 2.0 * t_prod(z, p.c, k))


Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p(order))
