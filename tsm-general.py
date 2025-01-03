#!/usr/bin/env python3
#  Example: ./general.py 15 4 0.010 10000 -1.0 0.0 1.0 -1.0 0.0 1.0 0.0 1.0 -1.0 1.0 -1.0 0.0 | ./plotAnimated.py -10 10
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_mul

def tsm_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])):
    pass

def get_p():
    return Parameters(a=float(argv[8]), b=float(argv[9]), c=float(argv[10]),
                      d=float(argv[11]), e=float(argv[12]), f=float(argv[13]),
                      g=float(argv[14]), h=float(argv[15]), i=float(argv[16]))

def ode(x, y, z, p, k):
    return Components(x=p.a * x[k] + p.b * y[k] + p.c * z[k],
                      y=p.d * x[k] + p.e * y[k] * p.f * z[k],
                      z=p.g * x[k] + p.h * y[k] + p.i * z[k])

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
