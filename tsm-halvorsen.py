#!/usr/bin/env python3
#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_sqr

class Parameters(namedtuple('ParametersType', ['a'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]))

def ode(x, y, z, p, k):
    return Components(x=-p.a * x[k] - 4.0 * (y[k] + z[k]) - t_sqr(y, k),
                      y=-p.a * y[k] - 4.0 * (z[k] + x[k]) - t_sqr(z, k),
                      z=-p.a * z[k] - 4.0 * (x[k] + y[k]) - t_sqr(x, k))

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
