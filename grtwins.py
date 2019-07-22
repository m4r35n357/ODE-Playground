#!/usr/bin/env python3

"""
Example: ./grtwins.py 2 -2 2 1001 0 1e-12 1e-12 | ./plotMany.py 2 300 >/dev/null

see https://www.mathpages.com/rr/s6-05/6-05.htm, equations (7) and (8)
"""

from sys import stderr, argv
from gmpy2 import get_context, sqrt, sin, cos, acos
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from playground import Mode, analyze
from taylor import to_mpfr


D1 = to_mpfr("1.0")
D15 = to_mpfr("1.5")
D2 = to_mpfr("2.0")
D3 = to_mpfr("3.0")

n = D1
r2 = to_mpfr("10.0")
q = r2 / D2
rtq = sqrt(q - D1)
pi = acos(to_mpfr("-1.0"))


def playground(a, value):
    return pi * n * (D1 + a.cos) ** D15 * q ** D15 - (
            q * rtq * (a + a.sin) + D2 * a * rtq + D2 * ((rtq + (a / D2).tan) / (rtq - (a / D2).tan)).ln) - value


def ratio(a):
    r1 = (D1 + cos(a)) * q
    return r1, ((a + sin(a)) * sqrt(q)) / (pi * n * (D1 + cos(a)) * sqrt(r1 - D3))


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

for result in analyze(playground, mode, x0, x1, steps, target, f_tol, x_tol):
    if result.mode == Mode.ROOT.name:
        print(result, file=stderr)
        answer = ratio(result.x)
        print(answer, file=stderr)
        print("n = {}".format(n), file=stderr)
        print("alpha = {}".format(result.x), file=stderr)
        print("r1 = {}".format(answer[0]), file=stderr)
        print("r2 = {}".format(r2), file=stderr)
        print("radial / circular = {}".format(answer[1]), file=stderr)
