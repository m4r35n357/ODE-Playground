#!/usr/bin/env python3

# Example:  ./solver.py NA 6 4.0 3.0 1e-12 1e-12 101 ROOT___

from sys import argv, exit, stderr
from ad import Context
from playground import Solver, Mode, bisect, newton

method = Solver.NA
try:
    method = Solver[argv[1]]
except KeyError:
    print(f"INVALID ANALYSIS: '{argv[1]}'", file=stderr)
    print(Solver.__members__.keys(), file=stderr)
    exit()
Context.places = int(argv[2])
xa = float(argv[3])
xb = float(argv[4])
f_tol = float(argv[5])
x_tol = float(argv[6])
max_it = int(argv[7])
mode = Mode.ROOT___
try:
    mode = Mode[argv[8]]
except KeyError:
    print(f"INVALID ANALYSIS MODE: '{argv[8]}'", file=stderr)
    print(Mode.__members__.keys(), file=stderr)
    exit()

xc = (xa + xb) / 2.0
model = lambda x: x**3 - 4 * x**2 + 3 * x - 2

if method in [Solver.BI, Solver.NA]:
    bisect(model, xa, xb, εf=f_tol, εx=x_tol, limit=max_it, mode=mode, debug=True)
    print("")
if method in [Solver.NT, Solver.NA]:
    newton(model, xc, εf=f_tol, εx=x_tol, limit=max_it, mode=mode, debug=True)
    print("")

# For Python console:
'''
from ad import *
from playground import *

model = lambda x: x**3 - 4 * x**2 + 3 * x - 2

x0 = 4.0
x1 = 3.0
x2 = (x0 + x1) / 2.0
y = 0.0

bisect(model, x0, x1, y=y, debug=True)
newton(model, x0, y=y, debug=True)

timeit(bisect(model, x0, x1, y=y))
timeit(newton(model, x2, y=y))

'''
