#!/usr/bin/env python3
#  Example: ./coupled.py 15 4 0.010 10000 -1.0 0.0 1.0 -1.0 0.0 1.0 | ./plotAnimated.py -10 10
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_mul

def tsm_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]), b=float(argv[9]), c=float(argv[10]))

def ode(x, y, z, p, k):
    return Components(x=p.b * y[k] + p.c * z[k],
                      y=p.c * z[k] * p.a * x[k],
                      z=p.a * x[k] + p.b * y[k])

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
