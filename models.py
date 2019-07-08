#!/usr/bin/env python3

from sys import argv
from playground import analyze
from functions import *

mode = int(argv[1])
assert mode == 0 or mode == 1 or mode == 2
x0 = float(argv[2])
x1 = float(argv[3])
assert x1 > x0
steps = int(argv[4])
assert steps > 0
order = int(argv[5])
assert order > 2
f_tol = float(argv[6])
x_tol = float(argv[7])

for result in analyze(composite1, mode, x0, x1, steps, f_tol, x_tol, max_it=1000, order=order):
    print(result, file=stderr)
