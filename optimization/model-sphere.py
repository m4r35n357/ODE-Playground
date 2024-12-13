#!/usr/bin/env python3
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
# CUT: ./model-sphere.py 8 256 100 0 10 >/dev/null
#  NM: ./model-sphere.py 8 10000 0 10 >/dev/null
from sys import argv
from random import random

def sphere (position):
    value = 0.0
    for j in range(len(position)):
        xo = position[j] - (j + 1)
        value += xo**2
    return value

if len(argv) == 6:
    from cut import coa
    dim = int(argv[2])
    m = int(argv[3])
    max_iter = int(argv[4])
    lower = float(argv[5])
    upper = float(argv[6])
    coa(sphere, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
    from nelder_mead import nelder_mead, Simplex
    dim = int(argv[2])
    tolerance = float(argv[3])
    max_iter = int(argv[4])
    lower = float(argv[5])
    upper = float(argv[6])
    nelder_mead(sphere, dim, max_iter, Simplex(sphere, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
    raise Exception('>>> Wrong number of arguments <<<')
