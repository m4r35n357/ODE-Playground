#!/usr/bin/env python3
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_const, t_mul

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c'])):
    pass

def get_p():
    return Parameters(a = float(argv[8]),
                      b = float(argv[9]),
                      c = float(argv[10]))

def ode(x, y, z, p, k):
    return Components(x = - y[k] - z[k],
                      y = x[k] + p.a * y[k],
                      z = t_const(p.b, k) + t_mul(x, z, k) - p.c * z[k])

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
