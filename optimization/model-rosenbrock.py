#!/usr/bin/env python3
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
# CUT: ./model-rosenbrock.py cut 8    256   100 0 10 >/dev/null
#  NM: ./model-rosenbrock.py  nm 8 1.0e-6 10000 0 10 >/dev/null
from sys import argv
from random import random

def rosenbrock (position):
    value = 0.0
    for j in range(len(position) - 1):
        xo = position[j] - (j + 1)
        xo_1 = position[j + 1] - (j + 2)
        value += (1.0 - xo)**2 + 100.0 * (xo_1 - xo**2)**2
    return value

if argv[1] == "cut":
    from cut import coa
    dim = int(argv[2])
    m = int(argv[3])
    max_iter = int(argv[4])
    lower = float(argv[5])
    upper = float(argv[6])
    coa(rosenbrock, dim, m, max_iter, lower, upper)
elif argv[1] == "nm":
    from nelder_mead import nelder_mead, Simplex
    dim = int(argv[2])
    tolerance = float(argv[3])
    max_iter = int(argv[4])
    lower = float(argv[5])
    upper = float(argv[6])
    nelder_mead(rosenbrock, dim, tolerance, max_iter, Simplex(rosenbrock, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
    raise Exception('>>> Wrong number of arguments <<<')
