#!/usr/bin/env python3

# Example:  ./solver.py NA -1.0  1.0 1e-15 1e-15 1000 2.0
# Example:  ./solver.py NA  1.0  3.0 1e-15 1e-15 1000 2.0

from sys import argv, exit, stderr
from playground import Solver, bisect, falsi, secant, newton, householder

method = Solver.NA
try:
    method = Solver[argv[1]]
except KeyError:
    print(f"INVALID ANALYSIS: '{argv[1]}'", file=stderr)
    print(Solver.__members__.keys(), file=stderr)
    exit()
xa = float(argv[2])
xb = float(argv[3])
f_tol = float(argv[4])
x_tol = float(argv[5])
max_it = int(argv[6])
target = float(argv[7])

model = lambda x: (x - 1)**2 / (x.cosh + 1).ln - 1
# model = lambda x: x**2
xc = (xa + xb) / 2.0
if (method == Solver.NA) or (method == Solver.BI):
    bisect(model, xa, xb, y=target, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.FP):
    falsi(model, xa, xb, y=target, εf=f_tol, εx=x_tol, limit=max_it, ill=False, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.FI):
    falsi(model, xa, xb, y=target, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.SC):
    secant(model, xa, xb, y=target, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.NT):
    newton(model, xc, y=target, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.H1):
    householder(model, xc, 2, y=target, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.H2):
    householder(model, xc, 3, y=target, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.H3):
    householder(model, xc, 4, y=target, εf=f_tol, εx=x_tol, limit=max, debug=True)
    print("")
if (method == Solver.NA) or (method == Solver.H4):
    householder(model, xc, 5, y=target, εf=f_tol, εx=x_tol, limit=max, debug=True)
    print("")


# For Python console:
'''
from ad import *
from playground import *

model = lambda x: (x - 1)**2 / (x.cosh + 1).ln - 1

x0 = 1.0
x1 = 3.0
x2 = (x0 + x1) / 2.0
y = 0.0

bisect(model, x0, x1, y=y, debug=True)
falsi(model, x0, x1, y=y, ill=False, debug=True)
falsi(model, x0, x1, y=y, debug=True)
secant(model, x0, x1, y=y, debug=True)
newton(model, x0, y=y, debug=True)
householder(model, x0, 2, y=y, debug=True)
householder(model, x0, 3, y=y, debug=True)
householder(model, x0, 4, y=y, debug=True)
householder(model, x0, 5, y=y, debug=True)

timeit(bisect(model, x0, x1, y=y))
timeit(falsi(model, x0, x1, y=y, ill=False))
timeit(falsi(model, x0, x1, y=y))
timeit(secant(model, x0, x1, y=y))
timeit(newton(model, x2, y=y))
timeit(householder(model, x2, 2, y=y))
timeit(householder(model, x2, 3, y=y))
timeit(householder(model, x2, 4, y=y))
timeit(householder(model, x2, 5, y=y))

'''
