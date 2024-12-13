#!/usr/bin/env python3
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
# CUT: ./model-treacle.py 8 256 100 0 10 >/dev/null
#  NM: ./model-treacle.py 8 10000 0 10 >/dev/null
from sys import argv
from math import sqrt
from random import random

def treacle (position):
    value = 0.0
    for j in range(len(position)):
        xo = position[j] - (j + 1)
        value += sqrt(abs(xo))
    return value

if len(argv) == 6:
    from cut import coa
    dim = int(argv[1])
    m = int(argv[2])
    max_iter = int(argv[3])
    lower = float(argv[4])
    upper = float(argv[5])
    coa(treacle, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
    from nelder_mead import nelder_mead, Simplex
    dim = int(argv[1])
    max_iter = int(argv[2])
    lower = float(argv[3])
    upper = float(argv[4])
    nelder_mead(treacle, dim, max_iter, Simplex(treacle, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
    raise Exception('>>> Wrong number of arguments <<<'
