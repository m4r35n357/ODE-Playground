#!/usr/bin/env python3
#  Example: ./zero.py 15 4 0.010 10000 -1.0 0.0 1.0 | ./plotAnimated.py -10 10
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from sys import argv, stderr
from collections import namedtuple
from ad import Components, Context, tsm, t_mul

def tsm_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

def get_p():
    return None

def ode(x, y, z, p, k):
    return Components(x=0.0,
                      y=0.0,
                      z=0.0)

Context.places, order, h, steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
tsm(ode, Context.places, order, h, steps, float(argv[5]), float(argv[6]), float(argv[7]), get_p())
