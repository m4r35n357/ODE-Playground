#!/usr/bin/env python3

from sys import stderr, argv
from gmpy2 import get_context
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from playground import analyze
from taylor import to_mpfr


def lorentz(a, value):
    # Example: ./models.py 0 .001 .999 1001 0 1e-12 1e-12 | ./plotMany.py 1 10 >/dev/null
    return 1 / (1 - a.sqr).sqrt - value


def cosx_x3(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return a.cos - a * a.sqr - value


def septic(a, value):
    # Example: ./models.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value


def composite1(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.exp + (a.sqr - 4.0).exp).ln - value


def composite2(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.sqr + (a.exp - 4).sqr).sqrt - value


def playground(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # return (2 * a).sin - 2 * a.sin * a.cos - value
    # return (2 * a).cos - a.cos.sqr + a.sin.sqr - value
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr - value
    # return a.sqr.sqrt - value
    # return a.exp.ln - value
    # return a.tan - a.sin / a.cos - value
    # return a.tanh - a.sinh / a.cosh - value
    # return a.sin.asin - value
    # return a.cos.acos - value
    return a.tan.atan - value


mode = int(argv[1])
assert mode == 0 or mode == 1 or mode == 2
x0 = to_mpfr(argv[2])
x1 = to_mpfr(argv[3])
assert x1 > x0
steps = int(argv[4])
assert steps > 0
target = to_mpfr(argv[5])
f_tol = to_mpfr(argv[6])
x_tol = to_mpfr(argv[7])

for result in analyze(composite1, mode, x0, x1, steps, target, f_tol, x_tol):
    print(result, file=stderr)
