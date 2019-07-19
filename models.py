#!/usr/bin/env python3

from sys import argv
from gmpy2 import get_context
get_context().precision = 236  # Set this BEFORE importing any AD stuff!
from ad import to_mpfr
from playground import analyze, Analysis
from functions import *

mode = int(argv[1])
assert mode in [item.value for item in Analysis.__members__.values()]
x0 = to_mpfr(argv[2])
x1 = to_mpfr(argv[3])
assert x1 > x0
steps = int(argv[4])
assert steps > 0
order = int(argv[5])
assert order > 2
f_tol = to_mpfr(argv[6])
x_tol = to_mpfr(argv[7])

for result in analyze(playground, mode, x0, x1, steps, f_tol, x_tol, max_it=1000, order=order):
    print(result, file=stderr)
