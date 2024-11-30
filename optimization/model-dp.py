#!/usr/bin/env python3
# CUT: ./model-dp.py cut 8    256   100 0 10 >/dev/null
#  NM: ./model-dp.py  nm 8 1.0e-6 10000 0 10 >/dev/null
from sys import argv
from random import random

def dixon_price (position):
    value = (position[0] - 1.0)**2
    for j in range(1, len(position)):
        value += (j + 1) * (2.0 * position[j]**2 - position[j - 1])**2
    return value

if len(argv) == 6:
    from cut import coa
    dim = int(argv[2])
    m = int(argv[3])
    max_iter = int(argv[4])
    lower = float(argv[5])
    upper = float(argv[6])
    coa(dixon_price, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
    from nelder_mead import nelder_mead, Simplex
    dim = int(argv[2])
    tolerance = float(argv[3])
    max_iter = int(argv[4])
    lower = float(argv[5])
    upper = float(argv[6])
    nelder_mead(dixon_price, dim, max_iter, Simplex(dixon_price, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
    raise Exception('>>> Wrong number of arguments <<<')
