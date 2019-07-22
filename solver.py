#!/usr/bin/env python3

# Example:  ./solver.py NA 2.0 0.0 1e-12 1000 2.0

from sys import argv, exit, stderr
from playground import Solver, bisect, secant, newton, householder

method = Solver.NA
try:
    method = Solver[argv[1]]
except KeyError:
    print('INVALID ANALYSIS', file=stderr)
    print(Solver.__members__.keys(), file=stderr)
    exit()
x0 = float(argv[2])
x1 = float(argv[3])
tol = float(argv[4])
max_it = int(argv[5])
target = float(argv[6])

# model = lambda x: x**2 / (x.cosh + 1).ln - 1
model = lambda x: x**2

if (method == Solver.NA) or (method == Solver.BI):
    print(bisect(model, x0, x1, εx=tol, εf=tol, limit=max_it, y=target, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.FP):
    print(secant(model, x0, x1, εx=tol, εf=tol, limit=max_it, y=target, fp=True, ill=False, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.FI):
    print(secant(model, x0, x1, εx=tol, εf=tol, limit=max_it, y=target, fp=True, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.SC):
    print(secant(model, x0, x1, εx=tol, εf=tol, limit=max_it, y=target, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.NT):
    print(newton(model, x0, εx=tol, εf=tol, limit=max_it, y=target, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.H1):
    print(householder(model, x0, 2, εx=tol, εf=tol, limit=max_it, y=target, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.H2):
    print(householder(model, x0, 3, εx=tol, εf=tol, limit=max_it, y=target, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.H3):
    print(householder(model, x0, 4, εx=tol, εf=tol, limit=max, y=target, debug=True))
    print("")
if (method == Solver.NA) or (method == Solver.H4):
    print(householder(model, x0, 5, εx=tol, εf=tol, limit=max, y=target, debug=True))
    print("")


# For Python console:
'''
from ad import *
from playground import *

model = lambda x: x**2 / (x.cosh + 1).ln - 1

x0 = 3.0
x1 = 0.0
y = 2.0

bisect(model, x0, x1, y=y, debug=True)
secant(model, x0, x1, y=y, fp=True, ill=False, debug=True)
secant(model, x0, x1, y=y, fp=True, debug=True)
secant(model, x0, x1, y=y, debug=True)
newton(model, x0, y=y, debug=True)
householder(model, x0, 2, y=y, debug=True)
householder(model, x0, 3, y=y, debug=True)
householder(model, x0, 4, y=y, debug=True)
householder(model, x0, 5, y=y, debug=True)

timeit(bisect(model, x0, x1, y=y))
timeit(secant(model, x0, x1, y=y, fp=True, ill=False))
timeit(secant(model, x0, x1, y=y, fp=True))
timeit(secant(model, x0, x1, y=y))
timeit(newton(model, x0, y=y))
timeit(householder(model, x0, 2, y=y))
timeit(householder(model, x0, 3, y=y))
timeit(householder(model, x0, 4, y=y))
timeit(householder(model, x0, 5, y=y))

'''
