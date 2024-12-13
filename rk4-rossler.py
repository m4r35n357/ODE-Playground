#!/usr/bin/env python3
#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, rk4

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]), b=float(argv[9]), c=float(argv[10]))

def ode(x, y, z, p):
    return Components(x = - y - z,
                      y = x + p.a * y,
                      z = p.b + (x - p.c) * z)

Context.places, skip, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
rk4(ode, Context.places, skip, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
