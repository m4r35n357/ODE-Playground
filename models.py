#!/usr/bin/env python3

from sys import argv, exit, stderr
from playground import analyze, Solver
from functions import playground

mode = Solver.NA
try:
    mode = Solver[argv[1]]
except KeyError:
    print('INVALID ANALYSIS', file=stderr)
    print(Solver.__members__.keys(), file=stderr)
    exit()
x0 = float(argv[2])
x1 = float(argv[3])
assert x1 > x0
steps = int(argv[4])
assert steps > 0
order = int(argv[5])
assert order > 2
f_tol = float(argv[6])
x_tol = float(argv[7])

for result in analyze(playground, mode, x0, x1, steps, f_tol, x_tol, limit=1000, order=order):
    print(result, file=stderr)
