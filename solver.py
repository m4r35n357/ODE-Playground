#!/usr/bin/env python3

from sys import argv, exit, stderr
from playground import Analysis, bisect, secant, newton, householder

method = Analysis.NA
try:
    method = Analysis[argv[1]]
except KeyError:
    print('INVALID ANALYSIS', file=stderr)
    print(Analysis.__members__.keys(), file=stderr)
    exit()
x0 = float(argv[2])
x1 = float(argv[3])
tol = float(argv[4])
max_it = int(argv[5])
target = float(argv[6])

model = lambda x: x**2

if (method == Analysis.NA) or (method == Analysis.BI):
    print(bisect(model, x0, x1, x_tol=tol, f_tol=tol, max_it=max_it, target=target, debug=True))
    print("")
if (method == Analysis.NA) or (method == Analysis.FP):
    print(secant(model, x0, x1, x_tol=tol, f_tol=tol, max_it=max_it, target=target, fp=True, debug=True))
    print("")
if (method == Analysis.NA) or (method == Analysis.SC):
    print(secant(model, x0, x1, x_tol=tol, f_tol=tol, max_it=max_it, target=target, debug=True))
    print("")
if (method == Analysis.NA) or (method == Analysis.NT):
    print(newton(model, x0, x_tol=tol, f_tol=tol, max_it=max_it, target=target, debug=True))
    print("")
if (method == Analysis.NA) or (method == Analysis.H1):
    print(householder(model, x0, 2, x_tol=tol, f_tol=tol, max_it=max_it, target=target, debug=True))
    print("")
if (method == Analysis.NA) or (method == Analysis.H2):
    print(householder(model, x0, 3, x_tol=tol, f_tol=tol, max_it=max_it, target=target, debug=True))
    print("")
if (method == Analysis.NA) or (method == Analysis.H3):
    print(householder(model, x0, 4, x_tol=tol, f_tol=tol, max_it=max, target=target, debug=True))
    print("")
